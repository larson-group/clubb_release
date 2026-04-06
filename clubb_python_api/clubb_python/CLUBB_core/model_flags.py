"""User-facing wrappers for routines from CLUBB_core/model_flags.F90."""

import clubb_f2py

from clubb_python.derived_types.config_flags import ConfigFlags, CONFIG_FLAG_NAMES
from clubb_python.derived_types.config_flags_converter import decode_config_flags, set_fortran_config_flags


def get_default_config_flags() -> ConfigFlags:
    """Return CLUBB default configuration flags as a ConfigFlags NamedTuple."""
    result = clubb_f2py.f2py_set_default_config_flags()
    return decode_config_flags(result[:len(CONFIG_FLAG_NAMES)])


def init_config_flags(flags: ConfigFlags):
    """Pack config flags into a Fortran-derived-type and store them in module storage."""
    set_fortran_config_flags(flags)
