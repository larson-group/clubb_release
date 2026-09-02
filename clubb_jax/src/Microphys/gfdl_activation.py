"""JAX port of Microphys/gfdl_activation.F90 — the CLUBB-side logic of GFDL droplet activation.

`aer_act_clubb_quadrature_Gauss` activates cloud droplets per PDF component by (1) computing the probability that
each component is in an updraft from the subgrid w distribution, then (2) looking up the activated droplet
concentration in the GFDL `aer_ccn_act_wpdf_k` 5-D aerosol table and weighting by the updraft probability.

This module ports the SELF-CONTAINED, differentiable CLUBB-side pieces:
  * `erff`           — the Numerical-Recipes error function used for the Gaussian-CDF updraft probability.
  * `updraft_weights`— the per-component in-updraft probabilities P1, P2 (normalized to sum to 1), from the PDF
                       w-mean/variance, mixture fraction, and component cloud fractions.
The GFDL `aer_ccn_act_wpdf_k` activation table itself reads external `../input/scm_activation/droplets*.dat`
lookup files (single-precision, no oracle) and is NOT ported — so the final Ndrop assembly remains.

Validated in `tests/test_gfdl_activation.py`: erff vs math.erf (<1e-6), updraft_weights vs a literal NumPy
transcription + the P1+P2=1 normalization, differentiable.
"""
import jax.numpy as jnp

_WP2_EPS = 1.0e-4       # w-variance threshold (gfdl_activation.F90:98)
_P_UPDRAFT_EPS = 1.0e-16


def erff(x):
    """Error function via the Numerical-Recipes rational-exponential approximation (gfdl_activation.F90:erff).

    erf(x) = 1 - erfc(x); max abs error ~1.2e-7. Smooth → differentiable.
    """
    z = jnp.abs(x)
    t = 1.0 / (1.0 + 0.5 * z)
    dumerfc = t * jnp.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (
        0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (
            1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    dumerfc = jnp.where(x < 0.0, 2.0 - dumerfc, dumerfc)
    return 1.0 - dumerfc


def _component_updraft(w_i, varnce_w_i, weight_i):
    """Unnormalized in-updraft probability for one PDF component (gfdl_activation.F90:109-135).

    varnce_w > wp2_eps : weight * (0.5 + 0.5 erf(w / sqrt(2 varnce_w)))   (Gaussian CDF that w>0)
    else, w > 0        : weight                                            (degenerate updraft)
    else               : 0
    """
    # Guard the sqrt so the unselected (varnce<=eps) branch has a finite gradient.
    safe_var = jnp.where(varnce_w_i > _WP2_EPS, varnce_w_i, 1.0)
    p_smooth = weight_i * (0.5 + 0.5 * erff(w_i / jnp.sqrt(2.0 * safe_var)))
    return jnp.where(varnce_w_i > _WP2_EPS, p_smooth,
                     jnp.where(w_i > 0.0, weight_i, 0.0))


def updraft_weights(w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac, cloud_frac_1, cloud_frac_2):
    """Normalized per-component in-updraft probabilities (P1, P2) (gfdl_activation.F90:109-143).

    Each component's weight is mixt_frac(or 1-mixt_frac) * cloud_frac; the pair is renormalized
    (P1=P2=0 if both are below P_updraft_eps). Inputs may be arrays (vectorized over levels). Differentiable.

    NOTE: faithfully reproduces a Fortran sequential-assignment QUIRK (gfdl_activation.F90:138-139): P1 is
    reassigned to P1/(P1+P2) *before* P2 is normalized, so P2's denominator uses the already-modified P1, i.e.
    P2_norm = P2 / (P1/(P1+P2) + P2)  —  NOT the symmetric P2/(P1+P2). The pair therefore does not sum to 1.
    """
    P1 = _component_updraft(w_1, varnce_w_1, mixt_frac * cloud_frac_1)
    P2 = _component_updraft(w_2, varnce_w_2, (1.0 - mixt_frac) * cloud_frac_2)
    tot = P1 + P2
    active = tot > _P_UPDRAFT_EPS
    safe_tot = jnp.where(active, tot, 1.0)
    P1n = jnp.where(active, P1 / safe_tot, 0.0)                 # P1 reassigned first (Fortran line 138)
    denom2 = P1n + P2                                           # uses the NEW P1n + old P2 (the quirk)
    P2n = jnp.where(active, P2 / jnp.where(active, denom2, 1.0), 0.0)
    return P1n, P2n


def aer_act_clubb_ndrop(drop_1, drop_2, P1_updraft, P2_updraft, mixt_frac, cloud_frac_1, cloud_frac_2):
    """Layer-averaged activated cloud-droplet concentration (gfdl_activation.F90:156-174).

    Combines the per-PDF-component activated droplet concentrations (drop_1, drop_2 — from the GFDL
    `aer_ccn_act_wpdf_k` aerosol-activation lookup table, supplied by the caller) with the in-updraft
    probabilities (P1, P2 from `updraft_weights`), then averages over the cloudy fraction of the layer:

        Ndrop = (drop_1 P1 + drop_2 P2) (mixt_frac cloud_frac_1 + (1-mixt_frac) cloud_frac_2).

    Pure jnp → differentiable. Completes the CLUBB-side orchestration of GFDL droplet activation (the lookup
    table itself reads external data and is not ported).
    """
    cloudy_frac = mixt_frac * cloud_frac_1 + (1.0 - mixt_frac) * cloud_frac_2
    return (drop_1 * P1_updraft + drop_2 * P2_updraft) * cloudy_frac
