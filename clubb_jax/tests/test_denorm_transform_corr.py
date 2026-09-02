#!/usr/bin/env python3
"""test_denorm_transform_corr.py — validate the JAX denorm_transform_corr port.

denorm_transform_corr (setup_clubb_pdf_params.F90:3208) transforms the normal-space PDF correlation arrays to
real ("standard") space: normal-normal pairs (chi/eta/w) carry over unchanged, normal-lognormal pairs use
corr_NN2NL, lognormal-lognormal pairs use corr_NN2LL. Oracle:
  1. Structural invariants (symmetric, unit diagonal).
  2. The three normal-normal entries equal the normal-space input exactly.
  3. Every transformed entry matches corr_NN2NL / corr_NN2LL called directly with the correct sigma/ratio args,
     including the Fortran quirk that component-2 Ncn transforms reuse the component-1 Ncn variance ratio.
  4. A finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import denorm_transform_corr
from clubb_jax.src.CLUBB_core.pdf_utilities import corr_NN2NL, corr_NN2LL

CHI, ETA, W, NCN = 0, 1, 2, 3
PDF_DIM = 6        # chi, eta, w, Ncn, rr, Nr
HM = [4, 5]
NG, NZT = 2, 4


def _inputs(seed):
    rng = np.random.default_rng(seed)
    # Symmetric unit-diagonal normal-space corr arrays.
    def _sym():
        L = rng.uniform(-0.7, 0.7, (NG, NZT, PDF_DIM, PDF_DIM))
        L = np.tril(L, -1)
        A = L + np.swapaxes(L, -1, -2)
        i = np.arange(PDF_DIM)
        A[:, :, i, i] = 1.0
        return A
    Cn1, Cn2 = _sym(), _sym()
    s1n = rng.uniform(0.1, 1.0, (NG, NZT, PDF_DIM))
    s2n = rng.uniform(0.1, 1.0, (NG, NZT, PDF_DIM))
    r1 = rng.uniform(0.05, 2.0, (NG, NZT, PDF_DIM))
    r2 = rng.uniform(0.05, 2.0, (NG, NZT, PDF_DIM))
    return s1n, s2n, r1, r2, Cn1, Cn2


def _call(s1n, s2n, r1, r2, Cn1, Cn2):
    return denorm_transform_corr(s1n, s2n, r1, r2, Cn1, Cn2, CHI, ETA, W, NCN)


def test_structure_and_normal_pairs():
    s1n, s2n, r1, r2, Cn1, Cn2 = _inputs(1)
    A1, A2 = _call(s1n, s2n, r1, r2, Cn1, Cn2)
    for A in (np.asarray(A1), np.asarray(A2)):
        assert np.max(np.abs(A - np.swapaxes(A, -1, -2))) < 1e-14, "not symmetric"
        assert np.max(np.abs(A[:, :, np.arange(PDF_DIM), np.arange(PDF_DIM)] - 1.0)) < 1e-14, "diag != 1"
    A1 = np.asarray(A1)
    for (rr, cc) in ((ETA, CHI), (W, CHI), (W, ETA)):
        assert np.max(np.abs(A1[:, :, rr, cc] - Cn1[:, :, rr, cc])) < 1e-14, "normal pair changed"
    print("  structure (symmetric, unit diagonal) + normal-normal pairs unchanged  PASS")


def test_transform_entries():
    s1n, s2n, r1, r2, Cn1, Cn2 = _inputs(2)
    A1, A2 = _call(s1n, s2n, r1, r2, Cn1, Cn2)
    A1 = np.asarray(A1); A2 = np.asarray(A2)

    # chi/eta/w x Ncn (NN2NL); comp 2 reuses r1[Ncn].
    for x in (CHI, ETA, W):
        e1 = corr_NN2NL(Cn1[:, :, NCN, x], s1n[:, :, NCN], r1[:, :, NCN])
        e2 = corr_NN2NL(Cn2[:, :, NCN, x], s2n[:, :, NCN], r1[:, :, NCN])
        assert np.max(np.abs(A1[:, :, NCN, x] - np.asarray(e1))) < 1e-13, f"Ncn-{x} comp1"
        assert np.max(np.abs(A2[:, :, NCN, x] - np.asarray(e2))) < 1e-13, f"Ncn-{x} comp2 (r1 quirk)"

    # chi/eta/w x hydromet (NN2NL).
    for x in (CHI, ETA, W):
        for j in HM:
            e1 = corr_NN2NL(Cn1[:, :, j, x], s1n[:, :, j], r1[:, :, j])
            e2 = corr_NN2NL(Cn2[:, :, j, x], s2n[:, :, j], r2[:, :, j])
            assert np.max(np.abs(A1[:, :, j, x] - np.asarray(e1))) < 1e-13, f"{j}-{x} comp1"
            assert np.max(np.abs(A2[:, :, j, x] - np.asarray(e2))) < 1e-13, f"{j}-{x} comp2"

    # Ncn x hydromet (NN2LL); comp 2 reuses r1[Ncn].
    for j in HM:
        e1 = corr_NN2LL(Cn1[:, :, j, NCN], s1n[:, :, NCN], s1n[:, :, j], r1[:, :, NCN], r1[:, :, j])
        e2 = corr_NN2LL(Cn2[:, :, j, NCN], s2n[:, :, NCN], s2n[:, :, j], r1[:, :, NCN], r2[:, :, j])
        assert np.max(np.abs(A1[:, :, j, NCN] - np.asarray(e1))) < 1e-13, f"Ncn-{j} comp1"
        assert np.max(np.abs(A2[:, :, j, NCN] - np.asarray(e2))) < 1e-13, f"Ncn-{j} comp2 (r1 quirk)"

    # hydromet x hydromet (NN2LL).
    e1 = corr_NN2LL(Cn1[:, :, 5, 4], s1n[:, :, 4], s1n[:, :, 5], r1[:, :, 4], r1[:, :, 5])
    e2 = corr_NN2LL(Cn2[:, :, 5, 4], s2n[:, :, 4], s2n[:, :, 5], r2[:, :, 4], r2[:, :, 5])
    assert np.max(np.abs(A1[:, :, 5, 4] - np.asarray(e1))) < 1e-13, "rr-Nr comp1"
    assert np.max(np.abs(A2[:, :, 5, 4] - np.asarray(e2))) < 1e-13, "rr-Nr comp2"
    print("  transform entries vs direct corr_NN2NL/corr_NN2LL (incl. Ncn-comp2 r1 quirk)  PASS")


def test_differentiable():
    s1n, s2n, r1, r2, Cn1, Cn2 = _inputs(3)
    def loss(C):
        A1, A2 = denorm_transform_corr(s1n, s2n, r1, r2, C, Cn2, CHI, ETA, W, NCN)
        return jnp.sum(A1 ** 2) + jnp.sum(A2 ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(Cn1)))
    assert np.isfinite(g).all(), "non-finite grad through denorm_transform_corr"
    print(f"  jax.grad through denorm_transform_corr: finite ({g.size} entries)  PASS")


def main():
    print("test_denorm_transform_corr:")
    for t in (test_structure_and_normal_pairs, test_transform_entries, test_differentiable):
        t()
    print("All denorm_transform_corr checks PASSED")


if __name__ == "__main__":
    main()
