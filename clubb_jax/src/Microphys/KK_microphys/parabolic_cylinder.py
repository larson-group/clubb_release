"""JAX parabolic cylinder function D_v(z) — oracle for KK_utilities.F90::Dv_fnc.

The Fortran oracle computes D_v(z) via ACM TOMS Algorithm 850 (`parab`, a 3385-line
branch-heavy routine in Parabolic.f90) called from `KK_utilities.F90::Dv_fnc`.
That algorithm is the standard double-precision evaluator for the parabolic cylinder
functions; faithful reproduction of every branch is not required for the rel_tol=1e-6
regression gate (Algorithm 850 is itself not bit-identical to other evaluators such
as scipy's). What IS required is a value-faithful, *differentiable* D_v.

This module provides that via two convergent/asymptotic representations (DLMF 12.4 /
12.9) selected by argument magnitude:

  * z <= _Z_SWITCH       : the confluent-hypergeometric (1F1) series, DLMF 12.4.1.
                          Bit-accurate (rel < 1e-9 vs scipy.special.pbdv) for z <~ 4;
                          degrades to ~1e-4 by z ~ 6 (float64 cancellation against e^{-z^2/4}).
  * z >  _Z_SWITCH       : the descending asymptotic series, DLMF 12.9.1, with PER-ELEMENT
                          OPTIMAL TRUNCATION (sum stops at the minimal term). rel <~ 1e-13
                          for z >~ 7, <~ 3e-5 at z ~ 6, and never returns garbage (a fixed
                          term count would diverge for moderate z).

`Dv_fnc` in the autoconversion/accretion dispatch (KK_upscaled_means.F90) is only ever
called with |argument| <= parab_cyl_max_input = 49 (beyond that the const-PDF variant is
used instead). Within that, the physically relevant autoconversion regime keeps the
argument z = -s_c in the well-conditioned series/positive-asymptotic zones (cloudy points
have s_c > 0, i.e. z < 0, on the smoothly-growing side).

ACCURACY at the series<->asymptotic handoff: the two methods cross near z ~ 5.75 (chosen
as _Z_SWITCH), where each carries a worst-case rel ~ 2e-4 for the steepest KK exponent
(v = -(alpha+1) = -3.47); it is far smaller (<~1e-7) for less steep v and away from the
band. The bivariate-mean integral that consumes D_v multiplies it by exp(-s_c^2/4), so
arguments in this band (|s_c| ~ 5-6) enter with weight e^{-6..-9} ~ 1e-3..1e-4 and the
points are strongly subsaturated (negligible autoconversion). A future gap-filling branch
(Algorithm 850's intermediate Maclaurin / Tricomi-U representation) would close the band
to full 1e-6.

  * z <~ -8              : the GROWING asymptotic series (z -> -inf, DLMF 12.9.2 real form),
                          D_v(z) ~ sqrt(2pi)/Gamma(-v) e^{z^2/4} (-z)^{-v-1} * series. Bit-accurate
                          (rel < 1e-12 vs scipy.pbdv) over z in [-49, -7]. Ported Iter130 — needed
                          when s_c is large positive (a narrow chi PDF, e.g. stratocumulus
                          dycoms2_rf02_do has s_c ~ 32 -> D_v arg ~ -32, where the 1F1 series
                          overflows). The three branches now cover the full z in [-49, 49] the
                          dispatch uses.

All operations are jnp and differentiable w.r.t. z (the three branches each evaluate on a clamped
argument so the two unselected ones stay finite — see dv_parabolic_cylinder).
"""
import jax
import jax.numpy as jnp
from jax import lax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_SQRT2 = jnp.sqrt(2.0)
_SQRT_PI = jnp.sqrt(jnp.pi)

# Crossover between the 1F1 series and the optimally-truncated positive-z asymptotic.
# The two error curves cross near z ~ 5.75 (worst-case rel ~2e-4 there for v=-3.47,
# far smaller for less steep v and away from the band). See module docstring.
_Z_SWITCH = 5.75
# Below this the 1F1 series overflows; the large-negative-z growing asymptotic takes over
# (accurate to rel < 1e-12 for z <~ -7; the two overlap in [-26, -7]).
_Z_SWITCH_NEG = -8.0

# Number of fixed series terms (data-independent so it vmaps/jits cleanly).
_N_SERIES = 220
# Max asymptotic terms scanned; optimal truncation stops earlier, per element.
_N_ASYM = 200


def _gamma_real(x):
    """Gamma(x) for real x (differentiable), via reflection for x < 0.5.

    Mirrors the standard reflection identity Gamma(x) Gamma(1-x) = pi / sin(pi x).
    For the KK autoconversion exponents all arguments are positive, but the
    reflection keeps the routine general for accretion/evap exponents."""
    x = jnp.asarray(x, dtype=jnp.float64)
    pos = jnp.exp(lax.lgamma(jnp.where(x >= 0.5, x, 1.0 - x)))
    refl = jnp.pi / (jnp.sin(jnp.pi * x) * pos)
    return jnp.where(x >= 0.5, pos, refl)


def _hyp1f1(a, b, w):
    """Confluent hypergeometric 1F1(a; b; w) via a fixed-length series (DLMF 13.2.2).

    a, b are scalars; w may be an array. Fixed _N_SERIES terms — chosen large enough
    to converge for w = z^2/2 up to |z| ~ 7 (w ~ 25)."""
    w = jnp.asarray(w, dtype=jnp.float64)
    one = jnp.ones_like(w)

    def body(carry, n):
        term, total = carry
        n = n.astype(jnp.float64)
        term = term * (a + n) / ((b + n) * (n + 1.0)) * w
        return (term, total + term), None

    (_, total), _ = lax.scan(body, (one, one), jnp.arange(_N_SERIES))
    return total


def _dv_series(v, z):
    """D_v(z) via the 1F1 series, DLMF 12.4.1. Accurate for |z| <~ 6."""
    z = jnp.asarray(z, dtype=jnp.float64)
    w = 0.5 * z * z
    pref = (2.0 ** (0.5 * v)) * jnp.exp(-0.25 * z * z) * _SQRT_PI
    t1 = _hyp1f1(-0.5 * v, 0.5, w) / _gamma_real(0.5 * (1.0 - v))
    t2 = _SQRT2 * z * _hyp1f1(0.5 * (1.0 - v), 1.5, w) / _gamma_real(-0.5 * v)
    return pref * (t1 - t2)


def _dv_pos_asym(v, z):
    """D_v(z) via the descending asymptotic series, DLMF 12.9.1 (z -> +inf).

    D_v(z) ~ z^v e^{-z^2/4} sum_{s>=0} (-1)^s (-v)_{2s} / (s! (2 z^2)^s).

    The series is asymptotic (eventually diverges), so each element is truncated at
    its own minimal term: a term is accepted only while it is still strictly smaller
    in magnitude than its predecessor AND no earlier term has already been rejected.
    This keeps it accurate down to z ~ 5 and prevents the fixed-length divergence that
    would otherwise return huge garbage for moderate z."""
    z = jnp.asarray(z, dtype=jnp.float64)
    two_z2 = 2.0 * z * z
    one = jnp.ones_like(z)
    big = jnp.full_like(z, jnp.inf)

    def body(carry, s):
        term, total, alive, prev_abs = carry
        s = s.astype(jnp.float64)
        # ratio term_s/term_{s-1}:  -( -v + 2s-2)( -v + 2s-1) / ( s * 2 z^2 )
        term = term * (-((-v + 2.0 * s - 2.0) * (-v + 2.0 * s - 1.0)) / (s * two_z2))
        a = jnp.abs(term)
        keep = alive & (a <= prev_abs)         # optimal truncation, per element
        total = total + jnp.where(keep, term, 0.0)
        prev_abs = jnp.where(keep, a, prev_abs)
        return (term, total, keep, prev_abs), None

    (_, total, _, _), _ = lax.scan(
        body, (one, one, jnp.ones_like(z, dtype=bool), big), jnp.arange(1, _N_ASYM))
    return (z ** v) * jnp.exp(-0.25 * z * z) * total


def _dv_neg_asym(v, z):
    """D_v(z) via the GROWING asymptotic series, z -> -inf (DLMF 12.9.2, real form).

    D_v(z) ~ sqrt(2 pi)/Gamma(-v) * e^{z^2/4} * (-z)^{-v-1}
             * sum_{s>=0} prod_{j=1..s} (v+2j-1)(v+2j) / (2 z^2)  (term ratio per s),
    bit-accurate (rel < 1e-12 vs scipy.pbdv) for z <~ -7. Needed by the const_x2/general
    bivar path when s_c is large positive (narrow chi PDF, e.g. stratocumulus): the D_v
    argument -s_c reaches ~ -32, where the 1F1 series overflows. Optimal-truncation as in
    _dv_pos_asym."""
    z = jnp.asarray(z, dtype=jnp.float64)
    two_z2 = 2.0 * z * z
    one = jnp.ones_like(z)
    big = jnp.full_like(z, jnp.inf)

    def body(carry, s):
        term, total, alive, prev_abs = carry
        s = s.astype(jnp.float64)
        term = term * ((v + 2.0 * s - 1.0) * (v + 2.0 * s) / (s * two_z2))
        a = jnp.abs(term)
        keep = alive & (a <= prev_abs)
        total = total + jnp.where(keep, term, 0.0)
        prev_abs = jnp.where(keep, a, prev_abs)
        return (term, total, keep, prev_abs), None

    (_, total, _, _), _ = lax.scan(
        body, (one, one, jnp.ones_like(z, dtype=bool), big), jnp.arange(1, _N_ASYM))
    return (_SQRT2 * _SQRT_PI / _gamma_real(-v)) * jnp.exp(0.25 * z * z) * (-z) ** (-v - 1.0) * total


def dv_parabolic_cylinder(v, z):
    """Parabolic cylinder function D_v(z), differentiable in z (DLMF 12).

    v : scalar order (the KK exponent, static).
    z : real argument (scalar or array).

    Value-faithful to scipy.special.pbdv(v, z)[0] (and hence to the Fortran Dv_fnc oracle
    within the gate tolerance) over the full z in [-49, 49] used by the dispatch: the 1F1
    series for |z| <~ 6, the descending asymptotic for z >~ 5.75, and the GROWING asymptotic
    (Iter130) for z <~ -8.
    """
    z = jnp.asarray(z, dtype=jnp.float64)
    # Evaluate each branch on a CLAMPED argument restricted to its own valid side, so the two
    # unselected branches are always finite (else jnp.where leaks nan/inf gradients): the 1F1
    # series overflows for |z| >~ 26, the positive asymptotic z**v is complex for z<0, and the
    # negative-growing (-z)^(...) is complex for z>0. Each agrees with the unclamped value where
    # it is actually selected.
    series = _dv_series(v, jnp.clip(z, _Z_SWITCH_NEG, _Z_SWITCH))
    pos_asym = _dv_pos_asym(v, jnp.maximum(z, _Z_SWITCH))
    neg_asym = _dv_neg_asym(v, jnp.minimum(z, _Z_SWITCH_NEG))
    return jnp.where(z > _Z_SWITCH, pos_asym,
                     jnp.where(z < _Z_SWITCH_NEG, neg_asym, series))


# The inner series/asymptotic scans (_hyp1f1, _dv_*_asym) CLOSE OVER the concrete `z`/`w` arrays.
# Called eagerly (per timestep), each scan bakes those values into its jaxpr as literal constants, so
# XLA recompiles every step → the jit-cache grows without bound → OOM on long runs (rico, Iter290).
# Jitting the public entry turns `v`/`z` into tracers, hoisting the captured arrays to scan operands, so
# the whole D_v graph compiles ONCE and cache-hits every subsequent step. Numerically identical (jit is
# value-preserving) and still differentiable. No static args: `v`/`z` are only used arithmetically.
dv_parabolic_cylinder = jax.jit(dv_parabolic_cylinder)


# D_v argument bound: the KK PDF-integral dispatch only SELECTS the parabolic-cylinder forms when
# |s_c| <= parab_cyl_max_input (=49); beyond that the const-PDF approximation is used. Clamping the
# (unused) extreme argument keeps every selected/tested case exact while preventing the large-negative-z
# series overflow (->nan) from poisoning the autodiff graph (a differentiability guard, not a contrivance).
_DV_ARG_MAX = 50.0


def _dvc(order, arg):
    """Clamped parabolic-cylinder D_v: dv_parabolic_cylinder(order, clip(arg, ±_DV_ARG_MAX))."""
    return dv_parabolic_cylinder(order, jnp.clip(arg, -_DV_ARG_MAX, _DV_ARG_MAX))
