"""JAX-side mirrors of `src/CLUBB_core` modules.

Porting deviations:
- This package marker has no Fortran counterpart. It exists only so Python can
  import the per-module JAX mirrors under `clubb_jax.src.CLUBB_core`.
"""

from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.err_info_codes import (
    ERR_NONE,
    ERR_XP2_XPYP_INVALID_C_UU,
    ERR_XP2_XPYP_MULTIPLE_LHS_REQUIRED,
)
from clubb_jax.src.CLUBB_core.grid_class import Grid
from clubb_jax.src.CLUBB_core.nu_vert_res_dep import NuVertResDep
from clubb_jax.src.CLUBB_core.pdf_params import implicit_coefs_terms, pdf_parameter
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx

__all__ = [
    "ErrInfo",
    "ERR_NONE",
    "ERR_XP2_XPYP_INVALID_C_UU",
    "ERR_XP2_XPYP_MULTIPLE_LHS_REQUIRED",
    "Grid",
    "NuVertResDep",
    "implicit_coefs_terms",
    "pdf_parameter",
    "SclrIdx",
]
