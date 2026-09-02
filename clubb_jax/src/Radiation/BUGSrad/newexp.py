"""Fast exp() approximation used by the BUGSrad two-stream solvers — port of newexp.F90.

BUGSrad's `module newexp` defines a `function exp(x)` that SHADOWS the intrinsic exp in the
two-stream radiative-transfer routines (two_rt_lw/sw and their _iter/_sel/_bs variants all
`use newexp`). It is a deliberately-cheap rational approximation, NOT the true exponential — so
to be BIT-FAITHFUL the JAX solvers must call this, NOT jnp.exp (cf. the KK epss=1e-4 parab and
the simplified-radiation low-accuracy paths: the Fortran's approximate special functions must be
replicated, not "improved").

  exp(x) ≈ 1 / (1 − x·(0.2507213 − x·(0.0292732 − x·0.0038278)))^4
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def newexp(x):
    """BUGSrad fast-exp approximation (newexp.F90). NOT the true exp — replicates the Fortran
    bit-for-bit so the two-stream solvers stay faithful."""
    x = jnp.asarray(x, dtype=jnp.float64)
    y1 = 1.0 - x * (0.2507213 - x * (0.0292732 - x * 0.0038278))
    y2 = y1 * y1
    y4 = y2 * y2
    return 1.0 / y4


# Exact Fortran name: `newexp.f90` declares this as `function exp` (it shadows the intrinsic exp
# inside `module newexp`). The JAX primary name is `newexp` to avoid shadowing jnp.exp at call
# sites; this module-scoped alias makes the bare Fortran name available (`from ...newexp import exp`,
# mirroring the solvers' `use newexp, only: exp`) — same convention as the jit-aliased raws.
exp = newexp
