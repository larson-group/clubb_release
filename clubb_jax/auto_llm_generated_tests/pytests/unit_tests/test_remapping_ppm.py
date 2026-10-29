#!/usr/bin/env python3
"""validate the JAX PPM (method 2) remapping port (remapping_module.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.remapping_module import remap_vals_ppm

_NG, _DZ, _ZTOP = 2, 50.0, 1500.0


def test_mass_conservation_refined():
    # Build a top→surface pressure grid; refine each cell into 2; PPM must conserve the pressure-weighted integral.
    rng = np.random.default_rng(8)
    ncol = 2
    # source edges surface→top (decreasing pressure)
    p_src = np.cumsum(np.concatenate([[1.0e5], -rng.uniform(2000, 4000, 20)]))[None, :].repeat(ncol, 0)
    # refined target: midpoints inserted
    mids = 0.5 * (p_src[:, :-1] + p_src[:, 1:])
    p_tgt = np.zeros((ncol, p_src.shape[1] + (p_src.shape[1] - 1)))
    p_tgt[:, ::2] = p_src
    p_tgt[:, 1::2] = mids
    src = rng.uniform(0.5, 3.0, (ncol, p_src.shape[1] - 1))
    tgt = np.asarray(remap_vals_ppm(p_src, p_tgt, src, iv=1))
    # pressure-weighted integral (mass) preserved: sum(src*dp_src) == sum(tgt*dp_tgt)
    dp_src = np.abs(np.diff(p_src, axis=1)); dp_tgt = np.abs(np.diff(p_tgt, axis=1))
    m_src = np.sum(src * dp_src, axis=1); m_tgt = np.sum(tgt * dp_tgt, axis=1)
    rel = np.max(np.abs(m_tgt - m_src) / np.abs(m_src))
    assert rel < 1e-12, f"PPM not mass-conservative on refined grid: rel {rel:.2e}"
    # finite + sane: kord=4 iv=1 PPM allows bounded edge overshoot (not strictly monotone), but the remapped
    # field must stay finite and within a generous band of the source range (no blow-up).
    assert np.isfinite(tgt).all(), "PPM produced non-finite values"
    span = src.max() - src.min()
    assert tgt.min() > src.min() - span and tgt.max() < src.max() + span, "PPM produced a runaway overshoot"
    print(f"  PPM mass conservation (refined grid): rel {rel:.2e}; finite + bounded  PASS")
