"""Type models and Fortran-derived-type storage converters for F2PY CLUBB APIs."""

from clubb_python.derived_types.common import Array
from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.config_flags import ConfigFlags, CONFIG_FLAG_NAMES
from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.hm_metadata import HmMetadata
from clubb_python.derived_types.pdf_params import (
    pdf_parameter,
    implicit_coefs_terms,
)

from clubb_python.derived_types.grid_class_converter import set_fortran_grid, get_fortran_grid
from clubb_python.derived_types.sclr_idx_converter import set_fortran_sclr_idx, get_fortran_sclr_idx
from clubb_python.derived_types.nu_vert_res_dep_converter import set_fortran_nu_vert_res_dep, get_fortran_nu_vert_res_dep
from clubb_python.derived_types.err_info_converter import init_err_info, set_fortran_err_info, get_fortran_err_info
from clubb_python.derived_types.pdf_params_converter import (
    get_fortran_pdf_params,
    set_fortran_pdf_params,
    get_fortran_pdf_params_zm,
    set_fortran_pdf_params_zm,
    get_fortran_implicit_coefs,
    set_fortran_implicit_coefs,
)
from clubb_python.derived_types.config_flags_converter import (
    decode_config_flags,
    encode_config_flags,
    set_fortran_config_flags,
)

__all__ = [
    'Array',
    'Grid',
    'ConfigFlags',
    'CONFIG_FLAG_NAMES',
    'SclrIdx',
    'NuVertResDep',
    'ErrInfo',
    'HmMetadata',
    'pdf_parameter',
    'implicit_coefs_terms',
    'set_fortran_grid',
    'get_fortran_grid',
    'set_fortran_sclr_idx',
    'get_fortran_sclr_idx',
    'set_fortran_nu_vert_res_dep',
    'get_fortran_nu_vert_res_dep',
    'init_err_info',
    'set_fortran_err_info',
    'get_fortran_err_info',
    'get_fortran_pdf_params',
    'set_fortran_pdf_params',
    'get_fortran_pdf_params_zm',
    'set_fortran_pdf_params_zm',
    'get_fortran_implicit_coefs',
    'set_fortran_implicit_coefs',
    'decode_config_flags',
    'encode_config_flags',
    'set_fortran_config_flags',
]
