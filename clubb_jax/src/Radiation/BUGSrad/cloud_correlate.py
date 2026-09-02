"""JAX port of Radiation/BUGSrad/cloud_correlate.F90 — generalized cloud-overlap (both subroutines).

OPTIONAL scheme — the default CLUBB/BUGSrad build (`-Dnooverlap`) does not use it — but part of the oracle.

`bugs_ctot` computes the column total cloud amount from a layer cloud-fraction profile and a vertical
decorrelation length `l_c`, interpolating between maximum overlap (l_c -> inf) and random overlap (l_c < 1 m).
The cloudy-layer condensation (overlap between consecutive CLOUDY layers, skipping clear gaps) is a discrete
boolean-mask selection; the overlap arithmetic and the `cld_below` recurrence — which the Fortran writes as a
nested triangular loop but is algebraically a CUMULATIVE PRODUCT — are jnp (`ctot_from_cloudy_layers`), so c_tot
is differentiable w.r.t. the cloud fractions within the multi-cloudy-layer branch.

`bugs_cloudfit` splits the profile into a maximally-overlapped part + a random part, finding the stacked cloud
fraction via a 1000-point grid search constrained by `c_tot`. That search + the threshold splits are DISCRETE
(non-differentiable) → ported as a faithful NumPy reference.

Validated in `tests/test_cloud_correlate.py`: both subroutines bit-faithful vs literal NumPy transcriptions of
the Fortran loops, the random-overlap analytic limit, the bugs_cloudfit reconstruction invariant
(c_maximal*cf_max + (1-c_maximal)*cf_random = acld), and an FD-correct grad of the bugs_ctot overlap core.
"""
import numpy as np
import jax.numpy as jnp

MIN_CF = 1.0e-6   # minimum layer cloud fraction considered cloudy
MIN_LC = 1.0      # minimum correlation length [m]; below this, random overlap is assumed
_R_DZ = 29.286    # R_d/g-ish constant in the hydrostatic dz estimate (cloud_correlate.F90:61)


def _midlayer_heights(pl2, tl):
    """Mid-layer heights from the hydrostatic relation (cloud_correlate.F90:59-64). Layers ordered 1=top..nlm=bot;
    z is built from the bottom up so z decreases with layer index."""
    nlm = tl.shape[0]
    dz = _R_DZ * jnp.log(pl2[1:nlm + 1] / pl2[0:nlm]) * tl   # dz per layer, shape (nlm,)
    # z0 accumulates from the bottom (index nlm-1) upward; z[k] = (running z0 below k) + dz[k]/2.
    dz_rev = dz[::-1]
    z0_below_rev = jnp.concatenate([jnp.zeros(1), jnp.cumsum(dz_rev)[:-1]])   # z0 before adding this layer
    z_rev = z0_below_rev + dz_rev / 2.0
    return z_rev[::-1]


def bugs_ctot(pl2, tl, acld, l_c):
    """Total cloud amount for a single column (cloud_correlate.F90:bugs_ctot, per i_domain).

    JAX is per-column (the Fortran `bugs_ctot` batches over the i_domain length); the leading
    domain-length args (nlen/len_loc/nlm) are dropped per the port's drop-batching-dims convention.

    Args:
        pl2: level pressure, shape (nlm+1,).
        tl:  layer temperature, shape (nlm,).
        acld: radiative cloud fraction, shape (nlm,).
        l_c: correlation length [m] (scalar).
    Returns:
        c_tot scalar.
    """
    acld = jnp.asarray(acld)
    cloudy = np.asarray(acld > MIN_CF)
    Ncloud = int(cloudy.sum())
    if Ncloud == 0:
        return jnp.asarray(0.0)
    if Ncloud == 1:
        return jnp.max(acld)
    if float(jnp.max(acld)) >= 1.0:
        return jnp.asarray(1.0)

    z = _midlayer_heights(jnp.asarray(pl2), jnp.asarray(tl))
    idx = np.nonzero(cloudy)[0]            # ascending cloudy-layer indices (DISCRETE selection)
    return ctot_from_cloudy_layers(acld[idx], z[idx], l_c)


def ctot_from_cloudy_layers(cloud, zc, l_c):
    """The differentiable core: total cloud amount from the condensed cloudy-layer fractions `cloud` and their
    heights `zc` (both shape (Ncloud,), Ncloud>=2) and correlation length `l_c`. Pure jnp → differentiable w.r.t.
    the cloud fractions (the discrete cloudy-layer SELECTION is done by the caller).

    cld_below[k] = (cloud[k]-olap[k-1]) * prod_{j<k-1} (1 - (cloud[j]-olap[j])/(1-cloud[j+1])); c_tot = sum.
    """
    cloud = jnp.asarray(cloud)
    Ncloud = cloud.shape[0]
    dz = zc[:-1] - zc[1:]
    if float(l_c) < MIN_LC:
        olap = cloud[:-1] * cloud[1:]      # pure random
    else:
        alpha = jnp.exp(-dz / l_c)
        olap = alpha * jnp.minimum(cloud[:-1], cloud[1:]) + (1.0 - alpha) * cloud[:-1] * cloud[1:]
    factor = 1.0 - (cloud[:-1] - olap) / (1.0 - cloud[1:])
    cp = jnp.cumprod(factor)
    cp_shift = jnp.concatenate([jnp.ones(1), cp])
    head = cloud[0:1]
    tail = (cloud[1:] - olap) * cp_shift[:Ncloud - 1]
    return jnp.sum(jnp.concatenate([head, tail]))


def bugs_cloudfit(acld, c_tot, tol=1.0e-5):
    """Split a column's cloud profile into a maximally-overlapped part + a random part
    (cloud_correlate.F90:bugs_cloudfit, per column; the Fortran's leading domain-length args are dropped).

    Given the total cloud amount `c_tot` (from bugs_ctot), find the stacked (maximally overlapped) cloud
    fraction `cf_stacked` such that the implied total matches `c_tot`, then split each layer's cloud fraction
    into `cf_max` (inside the maximal column) and `cf_random` (outside it), renormalized to [0,1].

    DISCRETE algorithm (a 1000-point grid search for cf_stacked + threshold splits) → NOT differentiable; ported
    as a faithful NumPy reference. Returns (c_maximal, cf_max, cf_random) as numpy arrays.

    Invariant (any nonzero column): c_maximal*cf_max[j] + (1-c_maximal)*cf_random[j] = acld[j].
    """
    acld = np.asarray(acld, dtype=np.float64)
    nlm = acld.shape[0]
    cloudy = acld > MIN_CF
    Nc = int(cloudy.sum())
    cf_max = np.zeros(nlm)
    cf_random = np.zeros(nlm)
    if Nc == 0:
        return 0.0, cf_max, cf_random
    min_frac = acld[cloudy].min()
    max_frac = acld[cloudy].max()

    c_tot_max = 1.0 - np.prod(1.0 - acld)   # purely random sky = maximum possible c_tot
    if c_tot_max - c_tot > tol:
        # Grid search: smallest cf_stacked = n*delta (n=0..1000) with c_tot_calcd(cf) <= c_tot.
        delta = max_frac / 1000.0
        n = np.arange(1001)
        cf = n * delta                                    # candidate cf_stacked values
        # c_tot_calcd(cf) = cf + (1-cf)*(1 - prod_{acld>cf}(1-acld)/(1-cf))
        above = acld[None, :] > cf[:, None]               # (1001, nlm)
        ratio = np.where(above, (1.0 - acld[None, :]) / (1.0 - cf[:, None]), 1.0)
        prod = np.prod(ratio, axis=1)
        c_tot_calcd = cf + (1.0 - cf) * (1.0 - prod)
        hit = np.nonzero(c_tot_calcd <= c_tot)[0]
        N = int(hit[0]) if hit.size else 1001             # first index meeting the exit condition
        # Fortran: loop sets cf_stacked=(N+1)*delta after exit; then if >max_frac -> max_frac, else -delta/2.
        cf_next = (N + 1) * delta
        cf_stacked = max_frac if cf_next > max_frac else cf_next - delta / 2.0
    else:
        # c_tot >= c_tot_max: minimum stacking. One cloud layer must stack at its own fraction.
        cf_stacked = min_frac if Nc == 1 else 0.0

    c_maximal = cf_stacked
    ge = acld >= cf_stacked
    cf_max = np.where(ge, cf_stacked, acld)
    cf_random = np.where(ge, acld - cf_stacked, 0.0)
    if 0.0 < c_maximal < 1.0:
        cf_max = cf_max / c_maximal
        cf_random = cf_random / (1.0 - c_maximal)
    return c_maximal, cf_max, cf_random
