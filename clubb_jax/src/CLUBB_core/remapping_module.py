"""JAX port of remapping_module.F90.

Porting deviations:
- The Fortran public ``remap_iv_from_var_name`` helper is not ported here; JAX
  callers currently pass ``iv`` directly.
- Debug-only conservation, monotonicity, consistency, and remap-matrix checks
  that only print/stop at debug level >= 2 are omitted from the jitted path.
- ``calc_mass_over_grid_intervals`` is vectorized over columns and integrates a
  cumulative piecewise-quadratic mass function rather than reproducing the
  Fortran sequential loop state.  It computes the same exact integral for the
  linear density spline, up to floating-point summation order.
- The Fortran optional ``R_ij`` reuse path is not exposed; JAX computes the
  remapping matrix functionally when the Ullrich method is selected.
- The PPM port implements the active ``kord = 4`` vertical remapping path.  The
  Fortran ``kord >= 7`` Huynh branch is not implemented because the CLUBB call
  site fixes ``kord`` to 4.
"""
from __future__ import annotations

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.CLUBB_core.constants_clubb import grav

_CONS_ULLRICH_REMAP = 1    # model_flags.F90 cons_ullrich_remap
_PPM_REMAP = 2             # model_flags.F90 ppm_remap

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#
#                GENERAL REMAPPING FUNCTIONS AND SUBROUTINES
#
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def _calc_mass_1col(rho_vals, rho_levels, grid_levels):
    """Single-column helper for ``calc_mass_over_grid_intervals``."""
    seg_mass = (rho_levels[1:] - rho_levels[:-1]) * (rho_vals[:-1] + rho_vals[1:]) * 0.5
    Mcum = jnp.concatenate([jnp.zeros((1,), dtype=rho_vals.dtype), jnp.cumsum(seg_mass)])  # mass node0→node j

    def M(z):
        # segment j with rho_levels[j] <= z <= rho_levels[j+1]
        j = jnp.clip(jnp.searchsorted(rho_levels, z, side='right') - 1, 0, rho_levels.shape[0] - 2)
        zl = rho_levels[j]; zh = rho_levels[j + 1]
        vl = rho_vals[j]; vh = rho_vals[j + 1]
        rho_z = vl + (vh - vl) * (z - zl) / (zh - zl)
        return Mcum[j] + (z - zl) * (vl + rho_z) * 0.5

    return M(grid_levels[1:]) - M(grid_levels[:-1])


def calc_mass_over_grid_intervals(rho_lin_spline_vals, rho_lin_spline_levels, grid_levels):
    """Calculate the mass over every interval (grid_levels(k) to grid_levels(k+1))
    of the grid levels (grid_levels).
    This function assumes a linear spline for rho, defined by the values of
    rho (lin_spline_rho_vals) given at the altitudes (lin_spline_rho_levels).
    So with the assumption that the linear spline is the exact rho function, this function
    computes the exact mass over each interval, such that any target grid with the same lowest
    and highest level would give the same total mass, if the linear spline is the same.

    Returns an array with the dimension grid_levels_idx-1 where at each index the mass of that
    interval is stored, starting at the bottom level.

    Args:
        rho_lin_spline_vals:   (ngrdcol, total_idx) density values [kg/m^3], ascending in altitude.
        rho_lin_spline_levels: (ngrdcol, total_idx) the corresponding altitudes [m], ascending.
        grid_levels:           (ngrdcol, nlev) grid level altitudes [m], ascending.

    Returns:
        mass_per_interval: (ngrdcol, nlev-1) mass over each cell [kg/m^2 per unit area].
    """
    rv = jnp.asarray(rho_lin_spline_vals, dtype=jnp.float64)
    rl = jnp.asarray(rho_lin_spline_levels, dtype=jnp.float64)
    gl = jnp.asarray(grid_levels, dtype=jnp.float64)
    # Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    #        lowest to highest in altitude
    #--------------------- Begin Code ---------------------
    # check if highest value of lin spline is higher or equal to highest level of target grid
    # find first level and set initial val_below and level_below
    # find all spline connection points between the level_below and the next grid level
    # in grid_levels and add up their masses, finally the rho value on the next grid level is
    # interpolated and the last part of mass added to the total mass of that interval
    return jax.vmap(_calc_mass_1col)(rv, rl, gl)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#                LINEAR REMAPPING METHOD BY ULLRICH ET AL.
#        formula (30) in 'Arbitrary-Order Conservative and Consistent
#                         Remapping and a Theory of Linear Maps: Part II'
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def remapping_matrix(levels_source, levels_target):
    """implements the remapping scheme proposed by Ullrich et al. in
    'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part II'
    in formula (30)

    Built from source/target *pressure* levels (eq. 30): for target cell k and source cell l, the overlap of
    the pressure intervals normalized by the target cell thickness:
      omega = |max(tgt[k+1],src[l+1]) − min(tgt[k],src[l])| if that max < that min, else 0;
      R[k,l] = omega / |tgt[k+1] − tgt[k]|.
    For identical source/target levels this is the identity (each cell overlaps only itself).

    Args:
        levels_source: (ngrdcol, ns) source pressure levels [Pa].
        levels_target: (ngrdcol, nt) target pressure levels [Pa].

    Returns:
        R_ij: (ngrdcol, nt-1, ns-1).
    """
    ls = jnp.asarray(levels_source, dtype=jnp.float64)
    lt = jnp.asarray(levels_target, dtype=jnp.float64)

    def _one(ls_i, lt_i):
        #--------------------- Begin Code ---------------------
        tgt_lo = lt_i[:-1][:, None]; tgt_hi = lt_i[1:][:, None]      # (nt-1, 1)
        src_lo = ls_i[:-1][None, :]; src_hi = ls_i[1:][None, :]      # (1, ns-1)
        omega_upper = jnp.maximum(tgt_hi, src_hi)
        omega_lower = jnp.minimum(tgt_lo, src_lo)
        omega = jnp.where(omega_upper < omega_lower, jnp.abs(omega_upper - omega_lower), 0.0)
        return omega / jnp.abs(tgt_hi - tgt_lo)

    return jax.vmap(_one)(ls, lt)


def matrix_vector_mult(x_vectors, A_matrices):
    """Calculates the matrix-vector product y = Ax, where A is a matrix with dim_output x dim_input
    and y,x are vectors with x having dim_input many entries and y having dim_output many entries.
    The result y is calculated for every ngrdcol

    Args:
        x_vectors:  (ngrdcol, dim_input).
        A_matrices: (ngrdcol, dim_output, dim_input).

    Returns:
        y_vectors: (ngrdcol, dim_output).
    """
    x = jnp.asarray(x_vectors, dtype=jnp.float64)
    A = jnp.asarray(A_matrices, dtype=jnp.float64)
    #--------------------- Begin Code ---------------------
    # matrix vector multiplication
    return jnp.einsum('ikj,ij->ik', A, x)


def _pressure_levels(mass_on_cells, p_sfc):
    """Hydrostatic pressure at each grid level: p[0]=p_sfc; p[k]=p[k-1] − mass[k-1]·grav. Vectorized over cols."""
    p_sfc = jnp.asarray(p_sfc, dtype=jnp.float64)
    # Calculate pressure levels from altitudes and densities
    dp = mass_on_cells * grav                                        # (ngrdcol, ncell)
    return p_sfc[:, None] - jnp.concatenate(
        [jnp.zeros((mass_on_cells.shape[0], 1), dtype=mass_on_cells.dtype), jnp.cumsum(dp, axis=1)], axis=1)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#                PIECEWISE PARABOLIC METHOD BY ULLRICH AND WOODWARD
#        in 'The Piecewise Parabolic Method (PPM) for gas-dynamical simulations'
#
# The implementation of the PPM was obtained from the Energy Exascale Earth
# System Model project, sponsored by the U.S.Department of Energy, Office of
# Science, Office of Biological and Environmental Research Earth Systems Model
# Development Program area of Earth and Environmental System Modeling.
#
# https://github.com/E3SM-Project/E3SM/blob/master/components/eam/src/dynamics/fv/mapz_module.F90
#
# The subroutines kmppm, steepz, ppm2m and map1_ppm were obtained
# from the E3SM code.
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def kmppm(dm, AA, AL, AR, A6, lmt):
    """the monotonicity constrained defined in the PPM paper."""
    if lmt == 3:
        return AL, AR, A6
    if lmt == 0:
        # Standard PPM constraint
        da1 = AR - AL
        da2 = da1 ** 2
        a6da = A6 * da1
        dm0 = dm == 0.0
        lo = (~dm0) & (a6da < -da2)
        hi = (~dm0) & (a6da > da2)
        A6n = jnp.where(dm0, 0.0, jnp.where(lo, 3.0 * (AL - AA), jnp.where(hi, 3.0 * (AR - AA), A6)))
        ARn = jnp.where(dm0, AA, jnp.where(lo, AL - 3.0 * (AL - AA), AR))   # AR = AL - A6 (lo branch)
        ALn = jnp.where(dm0, AA, jnp.where(hi, AR - 3.0 * (AR - AA), AL))   # AL = AR - A6 (hi branch)
        return ALn, ARn, A6n
    if lmt == 1:
        # Improved full monotonicity constraint (Lin)
        # Note: no need to provide first guess of A6 <-- a4(4,i)
        qmp = 2.0 * dm
        ALn = AA - jnp.sign(qmp) * jnp.minimum(jnp.abs(qmp), jnp.abs(AL - AA))
        ARn = AA + jnp.sign(qmp) * jnp.minimum(jnp.abs(qmp), jnp.abs(AR - AA))
        A6n = 3.0 * (2.0 * AA - (ALn + ARn))
        return ALn, ARn, A6n
    # Positive definite constraint
    cond = jnp.abs(AR - AL) < -A6
    a6_safe = jnp.where(A6 != 0.0, A6, 1.0)
    fmin = AA + 0.25 * (AR - AL) ** 2 / a6_safe + A6 / 12.0
    neg = cond & (fmin < 0.0)
    b1 = neg & (AA < AR) & (AA < AL)
    b2 = neg & (~b1) & (AR > AL)
    b3 = neg & (~b1) & (~b2)
    ALn = jnp.where(b1, AA, jnp.where(b3, AR - 3.0 * (AR - AA), AL))
    ARn = jnp.where(b1, AA, jnp.where(b2, AL - 3.0 * (AL - AA), AR))
    A6n = jnp.where(b1, 0.0, jnp.where(b2, 3.0 * (AL - AA), jnp.where(b3, 3.0 * (AR - AA), A6)))
    return ALn, ARn, A6n


def ppm2m(AA, delp, iv, kord):
    """PPM sub-grid reconstruction (remapping_module.F90:ppm2m), top→surface ordering. AA=cell means (ncol,km),
    delp=layer thicknesses (ncol,km). Returns (AL, AR, A6) for f(s)=AL+s[(AR-AL)+A6(1-s)]. Faithful to the
    kord=4 path (the only one map1_ppm uses); the kord>=7 Huynh branch is omitted."""
    ncol, km = AA.shape
    z = jnp.zeros((ncol, 1), dtype=AA.dtype)
    # Developer: S.-J. Lin, NASA-GSFC
    # Developer: A. Gettelman, Feb 2013, convert to column physics
    delq = AA[:, 1:] - AA[:, :-1]                      # (km-1): delq[j]=AA[j+1]-AA[j]
    d4 = delp[:, :-1] + delp[:, 1:]                    # (km-1): d4f[j]=delp[j]+delp[j+1] = Fortran d4(j+2)
    # Pad to align with Fortran d4(k)=delp(k-1)+delp(k): D4[k]=delp[k-1]+delp[k], k=1..km-1
    D4 = jnp.concatenate([z, d4], axis=1)              # (km): D4[k] valid for k>=1
    DQ = jnp.concatenate([delq, z], axis=1)            # (km): DQ[j]=AA[j+1]-AA[j] for j<=km-2

    # dc / df2 for interior k=1..km-2 (Fortran k=2..km1)
    AA_m1, AA_0, AA_p1 = AA[:, :-2], AA[:, 1:-1], AA[:, 2:]      # k-1, k, k+1
    c1 = (delp[:, :-2] + 0.5 * delp[:, 1:-1]) / D4[:, 2:]        # d4(k+1)
    c2 = (delp[:, 2:] + 0.5 * delp[:, 1:-1]) / D4[:, 1:-1]       # d4(k)
    tmp = delp[:, 1:-1] * (c1 * DQ[:, 1:-1] + c2 * DQ[:, :-2]) / (D4[:, 1:-1] + delp[:, 2:])
    qmax = jnp.maximum(jnp.maximum(AA_m1, AA_0), AA_p1) - AA_0
    qmin = AA_0 - jnp.minimum(jnp.minimum(AA_m1, AA_0), AA_p1)
    dc_int = jnp.sign(tmp) * jnp.minimum(jnp.minimum(jnp.abs(tmp), qmax), qmin)
    dc = jnp.concatenate([z, dc_int, z], axis=1)                 # (km), dc[0],dc[km-1] set later
    df2 = jnp.concatenate([z, tmp, z], axis=1)

    # 4th-order provisional left-edge AL for k=2..km-2 (Fortran k=3..km1)
    c1b = delq[:, 1:-1] * delp[:, 1:-2] / D4[:, 2:-1]            # delq(k-1)*delp(k-1)/d4(k)
    a1 = D4[:, 1:-2] / (D4[:, 2:-1] + delp[:, 1:-2])            # d4(k-1)/(d4(k)+delp(k-1))
    a2 = D4[:, 3:] / (D4[:, 2:-1] + delp[:, 2:-1])             # d4(k+1)/(d4(k)+delp(k))
    AL_int = AA[:, 1:-2] + c1b + 2.0 / (D4[:, 1:-2] + D4[:, 3:]) * (
        delp[:, 2:-1] * (c1b * (a1 - a2) + a2 * dc[:, 1:-2]) - delp[:, 1:-2] * a1 * dc[:, 2:-1])
    AL = jnp.zeros_like(AA)
    AL = AL.at[:, 2:-1].set(AL_int)

    # steepz discontinuity steepening (modifies AL interior)
    AL = steepz(AA, AL, df2, dc, DQ, delp, D4)

    AR = jnp.zeros_like(AA)
    # Top boundary cubic (k=1,2 Fortran → 0,1)
    d1 = delp[:, 0]; d2 = delp[:, 1]
    qm = (d2 * AA[:, 0] + d1 * AA[:, 1]) / (d1 + d2)
    dq = 2.0 * (AA[:, 1] - AA[:, 0]) / (d1 + d2)
    c1t = 4.0 * (AL[:, 2] - qm - d2 * dq) / (d2 * (2.0 * d2 * d2 + d1 * (d2 + 3.0 * d1)))
    c3t = dq - 0.5 * c1t * (d2 * (5.0 * d1 + d2) - 3.0 * d1 ** 2)
    AL1 = qm - 0.25 * c1t * d1 * d2 * (d2 + 3.0 * d1)          # AL(2)
    AL0 = d1 * (2.0 * c1t * d1 ** 2 - c3t) + AL1               # AL(1)
    dc0 = AA[:, 0] - AL0
    cmax = jnp.maximum(AA[:, 0], AA[:, 1]); cmin = jnp.minimum(AA[:, 0], AA[:, 1])
    AL1 = jnp.minimum(jnp.maximum(AL1, cmin), cmax)
    AL = AL.at[:, 0].set(AL0).at[:, 1].set(AL1)
    dc = dc.at[:, 0].set(dc0)
    if iv == 0:
        AL = AL.at[:, 0].set(jnp.maximum(0.0, AL[:, 0])).at[:, 1].set(jnp.maximum(0.0, AL[:, 1]))
    elif iv == -1:
        a10 = AA[:, 0]
        AL0w = jnp.where(a10 * AL[:, 0] <= 0.0, 0.0,
                         jnp.sign(a10) * jnp.minimum(jnp.abs(a10), jnp.abs(AL[:, 0])))
        AL = AL.at[:, 0].set(AL0w)

    # Bottom boundary cubic (k=km,km1 Fortran → km-1,km-2)
    d1b = delp[:, -1]; d2b = delp[:, -2]
    qmb = (d2b * AA[:, -1] + d1b * AA[:, -2]) / (d1b + d2b)
    dqb = 2.0 * (AA[:, -2] - AA[:, -1]) / (d1b + d2b)
    c1bm = (AL[:, -2] - qmb - d2b * dqb) / (d2b * (2.0 * d2b * d2b + d1b * (d2b + 3.0 * d1b)))
    c3bm = dqb - 2.0 * c1bm * (d2b * (5.0 * d1b + d2b) - 3.0 * d1b ** 2)
    ALkm = qmb - c1bm * d1b * d2b * (d2b + 3.0 * d1b)          # AL(km)
    ARkm = d1b * (8.0 * c1bm * d1b ** 2 - c3bm) + ALkm         # AR(km)
    dckm = ARkm - AA[:, -1]
    cmaxb = jnp.maximum(AA[:, -1], AA[:, -2]); cminb = jnp.minimum(AA[:, -1], AA[:, -2])
    ALkm = jnp.minimum(jnp.maximum(ALkm, cminb), cmaxb)
    AL = AL.at[:, -1].set(ALkm)
    AR = AR.at[:, -1].set(ARkm)
    dc = dc.at[:, -1].set(dckm)
    if iv == 0:
        AR = AR.at[:, -1].set(jnp.maximum(0.0, AR[:, -1]))
    elif iv == -1:
        a1km = AA[:, -1]
        ARkmw = jnp.where(a1km * AR[:, -1] <= 0.0, 0.0,
                          jnp.sign(a1km) * jnp.minimum(jnp.abs(a1km), jnp.abs(AR[:, -1])))
        AR = AR.at[:, -1].set(ARkmw)

    # AR(k) = AL(k+1) for k=1..km-1 (0-based 0..km-2); AR(km) already set above
    AR = AR.at[:, :-1].set(AL[:, 1:])

    # A6 + kmppm per level group. lmt=0 for the top-2 and bottom-2 layers, lmt=kord-3 for the interior.
    A6 = 3.0 * (2.0 * AA - (AL + AR))
    lmt_int = max(0, kord - 3)
    if iv == 0:
        lmt_int = min(2, lmt_int)
    # interior uses kmppm(lmt_int) which (for lmt 1/2) recomputes AL/AR/A6 internally.
    AL0g, AR0g, A60g = kmppm(dc, AA, AL, AR, A6, 0)
    ALig, ARig, A6ig = kmppm(dc, AA, AL, AR, A6, lmt_int)
    k = jnp.arange(km)
    interior = (k >= 2) & (k <= km - 3)                          # Fortran k=3..km-2
    sel = interior[None, :]
    AL = jnp.where(sel, ALig, AL0g)
    AR = jnp.where(sel, ARig, AR0g)
    A6 = jnp.where(sel, A6ig, A60g)
    return AL, AR, A6


def steepz(AA, AL, df2, dm, dq, dp, d4):
    """the discontinuity adjustment defined in the PPM paper."""
    ncol, km = AA.shape
    z = jnp.zeros((ncol, 1), dtype=AA.dtype)
    #--------------------- Begin Code ---------------------
    # Compute ratio of dq/dp
    # rat(k) = dq(k-1)/d4(k) for k=2..km (Fortran) → 0-based k=1..km-1: rat[k]=DQ[k-1]/D4[k]
    d4_safe = jnp.where(d4 != 0.0, d4, 1.0)
    rat = jnp.concatenate([z, dq[:, :-1] / d4_safe[:, 1:]], axis=1)   # rat[0]=0, rat[k]=DQ[k-1]/D4[k]
    # f(k) = (rat(k+1)-rat(k))/(dp(k-1)+dp(k)+dp(k+1)) for k=2..km-1 (0-based 1..km-2)
    f = jnp.zeros_like(AA)
    denom_f = dp[:, :-2] + dp[:, 1:-1] + dp[:, 2:]
    f = f.at[:, 1:-1].set((rat[:, 2:] - rat[:, 1:-1]) / denom_f)
    # alfa(k) for Fortran k=3..km-2 (0-based k=2..km-3, km-4 elements). f(k±1) etc. shift accordingly.
    fkm1, fk, fkp1 = f[:, 1:-3], f[:, 2:-2], f[:, 3:-1]         # f(k-1), f(k), f(k+1)
    df2k = df2[:, 2:-2]
    dg2 = (fkp1 - fkm1) * ((dp[:, 3:-1] - dp[:, 1:-3]) ** 2 + d4[:, 2:-2] * d4[:, 3:-1])
    df2_safe = jnp.where(df2k != 0.0, df2k, 1.0)
    alfa_int = jnp.where((fkp1 * fkm1 < 0.0) & (df2k != 0.0),
                         jnp.maximum(0.0, jnp.minimum(0.5, -3.0 / 16.0 * dg2 / df2_safe)), 0.0)
    alfa = jnp.zeros_like(AA)
    alfa = alfa.at[:, 2:-2].set(alfa_int)
    # AL(k) update for k=4..km-2 (0-based 3..km-3)
    a_km1, a_k = alfa[:, 2:-3], alfa[:, 3:-2]                   # alfa(k-1), alfa(k) for 0-based k=3..km-3
    AL_new = ((1.0 - a_km1 - a_k) * AL[:, 3:-2]
              + a_km1 * (AA[:, 3:-2] - dm[:, 3:-2])
              + a_k * (AA[:, 2:-3] + dm[:, 2:-3]))
    return AL.at[:, 3:-2].set(AL_new)


def map1_ppm(pe1, q1, pe2, iv, kord):
    """E3SM PPM vertical remap (remapping_module.F90:map1_ppm), all arrays TOP→SURFACE (increasing pressure).
    pe1 (ncol, km+1) source edges, q1 (ncol, km) source means, pe2 (ncol, kn+1) target edges. Returns q2
    (ncol, kn). Conservative: integrates the PPM sub-grid parabola of each source cell over each target cell."""
    ncol, kmp1 = pe1.shape
    km = kmp1 - 1
    kn = pe2.shape[1] - 1
    dp1 = pe1[:, 1:] - pe1[:, :-1]                              # (km)
    AL, AR, A6 = ppm2m(q1, dp1, iv, kord)
    AA = q1

    # k0(k): smallest source cell kk with pe2(k) <= pe1(kk+1). pe1 increasing → searchsorted.
    # For each target edge pe2[:,k], k0 = clip(searchsorted(pe1_row, pe2val, 'left')-1, 0, km-1).
    def _k0_row(pe1r, pe2r):
        idx = jnp.searchsorted(pe1r, pe2r, side='left') - 1
        return jnp.clip(idx, 0, km - 1)
    k0 = jax.vmap(_k0_row)(pe1, pe2)                            # (ncol, kn+1)

    # cumulative cell masses for the whole-cell contribution sum
    dpq = dp1 * AA                                              # (km)
    cum = jnp.concatenate([jnp.zeros((ncol, 1), dtype=AA.dtype), jnp.cumsum(dpq, axis=1)], axis=1)  # (km+1)

    def _gatherrow(M, idx):  # M (km..), idx (kn..) → gather per row done via vmap outside
        return M[idx]

    r3 = 1.0 / 3.0
    r23 = 2.0 / 3.0

    def _remap_row(pe1r, pe2r, dp1r, ALr, ARr, A6r, AAr, k0r, cumr):
        kk = k0r[:-1]                                           # (kn) source cell of pe2(k)
        kk2 = k0r[1:]                                           # (kn) source cell of pe2(k+1)
        dpkk = dp1r[kk]
        # Consider contribution between pe1(i,kk) and pe2(i,k)
        pl = (pe2r[:-1] - pe1r[kk]) / dpkk
        # Check to see if pe2(i,k+1) and pe2(i,k) are in same pe1 interval
        pr = (pe2r[1:] - pe1r[kk]) / dpkk
        q_same = ALr[kk] + 0.5 * (A6r[kk] + ARr[kk] - ALr[kk]) * (pr + pl) - A6r[kk] * r3 * (pr * (pr + pl) + pl ** 2)
        # Consider contribution between pe2(i,k) and pe1(i,kk+1)
        qsum1 = (pe1r[kk + 1] - pe2r[:-1]) * (
            ALr[kk] + 0.5 * (A6r[kk] + ARr[kk] - ALr[kk]) * (1.0 + pl) - A6r[kk] * (r3 * (1.0 + pl * (1.0 + pl))))
        # Next consider contribution between pe1(i,kk+1) and pe1(i,k0(i,k+1))
        qsum_mid = cumr[jnp.maximum(kk2, 0)] - cumr[jnp.minimum(kk + 1, km)]   # sum dp1*AA over (kk+1 .. kk2-1)
        klb = kk2
        delp = pe2r[1:] - pe1r[klb]
        esl = delp / dp1r[klb]
        # Now consider contribution between pe1(i,k0(i,k+1)) and pe2(i,k+1)
        qsum3 = delp * (ALr[klb] + 0.5 * esl * (ARr[klb] - ALr[klb] + A6r[klb] * (1.0 - r23 * esl)))
        q_multi = (qsum1 + qsum_mid + qsum3) / (pe2r[1:] - pe2r[:-1])
        return jnp.where(kk2 == kk, q_same, q_multi)

    return jax.vmap(_remap_row)(pe1, pe2, dp1, AL, AR, A6, AA, k0, cum)


def remap_vals_ppm(levels_source, levels_target, source_values, iv, kord=4):
    """This is a wrapper function for E3SM's subroutine implementing the
    Piecewise Parabolic Method (PPM) that takes in the pressure levels from surface to top
    """
    ls = jnp.asarray(levels_source, dtype=jnp.float64)
    lt = jnp.asarray(levels_target, dtype=jnp.float64)
    sv = jnp.asarray(source_values, dtype=jnp.float64)
    # the pressure levels and values are flipped, since map1_ppm takes in the grids
    # from top to surface instead of surface to top, like CLUBB does
    pe1 = ls[:, ::-1]                                           # top→surface (increasing pressure)
    pe2 = lt[:, ::-1]
    q1 = sv[:, ::-1]
    # remap values using E3SM's PPM implementation
    q2 = map1_ppm(pe1, q1, pe2, iv, kord)
    # flip output values back so they correspond to the surface to top ordering in CLUBB's grid
    return q2[:, ::-1]                                          # back to surface→top


def remap_vals_to_target(levels_source, levels_target, rho_lin_spline_vals, rho_lin_spline_levels,
                         source_values, p_sfc, grid_remap_method=_CONS_ULLRICH_REMAP, iv=1):
    """Remap source_values (defined on the source-grid cells) to the target-grid cells, conserving mass
    (remapping_module.F90:remap_vals_to_target_helper, Ullrich-linear branch).

    Args:
        levels_source: (ngrdcol, ns) source grid level altitudes [m], ascending.
        levels_target: (ngrdcol, nt) target grid level altitudes [m], ascending.
        rho_lin_spline_vals/levels: (ngrdcol, total_idx) density spline.
        source_values: (ngrdcol, ns-1) field on source cells.
        p_sfc: (ngrdcol,) surface pressure [Pa].
        grid_remap_method: 1 = Ullrich linear (only method implemented).

    Returns:
        target_values: (ngrdcol, nt-1).
    """
    if grid_remap_method not in (_CONS_ULLRICH_REMAP, _PPM_REMAP):
        raise ValueError(f"no remapping method for grid_remap_method={grid_remap_method}")

    # Calculate pressure levels from altitudes and densities
    mass_src = calc_mass_over_grid_intervals(rho_lin_spline_vals, rho_lin_spline_levels, levels_source)
    mass_tgt = calc_mass_over_grid_intervals(rho_lin_spline_vals, rho_lin_spline_levels, levels_target)
    p_src = _pressure_levels(mass_src, p_sfc)
    p_tgt = _pressure_levels(mass_tgt, p_sfc)
    # Remap values with the selected remapping method
    if grid_remap_method == _PPM_REMAP:
        return remap_vals_ppm(p_src, p_tgt, source_values, iv)
    # Construct remapping matrix for cons_ullrich_remap method, if not
    # already given as an optional parameter
    R = remapping_matrix(p_src, p_tgt)
    return matrix_vector_mult(source_values, R)


def remap_vals_to_target_same_grid(source_values, zm, zt, rho_lin_spline_vals, rho_lin_spline_levels,
                                   p_sfc, grid_remap_method=_CONS_ULLRICH_REMAP, l_zt_variable=True, iv=1):
    """f2py-exposed driver: remap within a single (source == target) grid (remapping_module.F90:
    remap_vals_to_target, with gr_source == gr_target == stored grid). For a same grid the remapping matrix
    is the identity, so the output equals source_values — but the full mass/pressure/matrix pipeline is run
    faithfully.

    Args:
        source_values: (ngrdcol, ncell) field on the source cells (zt cells if l_zt_variable else zm cells).
        zm: (ngrdcol, nzm) momentum-level altitudes [m].
        zt: (ngrdcol, nzt) thermodynamic-level altitudes [m]  (nzt = nzm-1).
        rho_lin_spline_vals/levels: (ngrdcol, total_idx) density spline.
        p_sfc: (ngrdcol,) surface pressure [Pa].
        grid_remap_method: 1 = Ullrich linear.
        l_zt_variable: True if source_values live on zt levels (use zm as cell edges); False if on zm levels
                       (cell edges = [zm[0], zt..., zm[-1]]).

    Returns:
        target_values: (ngrdcol, ncell), same shape as source_values.
    """
    zm = jnp.asarray(zm, dtype=jnp.float64)
    zt = jnp.asarray(zt, dtype=jnp.float64)
    if l_zt_variable:
        # if variable is given on the zt levels, just take the surrounding zm levels to build
        # the corresponding grid cell
        levels = zm
    else:
        # construct surrounding grid levels for values given on the zm levels
        # for inner levels take the surrounding zt levels, for the outer levels,
        # take the zt level and the zm level
        levels = jnp.concatenate([zm[:, :1], zt, zm[:, -1:]], axis=1)
    return remap_vals_to_target(levels, levels, rho_lin_spline_vals, rho_lin_spline_levels,
                                source_values, p_sfc, grid_remap_method, iv=iv)
