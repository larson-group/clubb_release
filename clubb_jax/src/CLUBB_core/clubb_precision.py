"""Floating-point precision of CLUBB-JAX, the counterpart of clubb_precision.F90.

CLUBB-JAX runs either in double precision (the default) or in single
precision. The Fortran build fixes its working real kind, core_rknd, at
compile time. JAX instead uses one process-wide setting, jax_enable_x64: when
it is on, arrays created without an explicit dtype are float64; when it is off,
every JAX array is float32. CLUBB-JAX uses that setting as its core_rknd, so
model code creates arrays without a dtype and follows the selected precision.

Usage:

* Set CLUBB_JAX_PRECISION before starting Python, for example
  ``CLUBB_JAX_PRECISION=single ./run_scripts/run_scm.py -jax bomex``. Unset
  means double; ``clubb_jax.backends.common.precision_name`` lists the accepted
  spellings. The launcher uses the same function to report the selection
  without importing JAX.
* ``configure_jax_precision()`` applies the selection to JAX. The
  ``clubb_jax.src`` package calls it once on import, before any submodule
  creates arrays, so model code and tests do not call it. Import
  ``clubb_jax.src`` before creating JAX arrays in a script, and do not set
  jax_enable_x64 yourself: arrays created earlier keep their dtype, and the
  package import overrides an earlier setting.
* ``core_rknd`` is the JAX dtype of the working precision. Use it for
  precision queries, such as ``jnp.finfo(core_rknd).eps`` for Fortran
  ``epsilon(1._core_rknd)``, and to convert integer or boolean values, such
  as ``mask.astype(core_rknd)`` for Fortran ``real(x, kind=core_rknd)``.
  Arrays created from float values need no dtype. It is chosen from
  USE_X64 rather than fixed to ``jnp.float64``: JAX turns a requested
  float64 into float32 when x64 is off, but ``jnp.finfo(jnp.float64)``
  still reports double precision.
* ``CLUBB_JAX_PRECISION`` is the normalized selection, ``"single"`` or
  ``"double"``, and ``USE_X64`` is True for double precision.

Fixed precision:

The Fortran also uses kinds that do not follow core_rknd. CLUBB-JAX keeps two
of them because their results feed the JAX-versus-Fortran comparisons:

* Morrison microphysics computes in float32, as its default-real Fortran does.
* The inverse-normal coefficients in SILHS ``cdfnorminv`` are rounded to
  float32, as the Fortran default-real literals are.

Model times (time_initial, dt_main, time_current, ...) are Python floats,
which are double precision in both modes, like Fortran real(time_precision).

Single precision:

Single precision is for performance and memory experiments, mainly on GPUs.
Apart from the fixed-precision code above, all JAX arrays are float32,
including code the Fortran keeps in double precision (``dp``) in any build,
such as pdf_utilities, the KK microphysics integrals, BUGSrad radiation and
statistics accumulation. Its results therefore do not match a single-precision
Fortran build bit for bit. On bomex, single-precision CLUBB-JAX and Fortran
differ from their double-precision runs by similar amounts, but
dycoms2_rf02_do (KK microphysics) produces NaNs in single precision.
"""

import os

import jax
import jax.numpy as jnp

from clubb_jax.backends.common import precision_name


CLUBB_JAX_PRECISION = precision_name(os.environ.get("CLUBB_JAX_PRECISION"))
USE_X64 = CLUBB_JAX_PRECISION == "double"

# Working real kind; see "Usage" above.
core_rknd = jnp.float64 if USE_X64 else jnp.float32


def configure_jax_precision() -> None:
    """Apply the process-wide precision selected by CLUBB_JAX_PRECISION."""
    jax.config.update("jax_enable_x64", USE_X64)
