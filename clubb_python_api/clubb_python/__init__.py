"""User-facing Python package for CLUBB F2PY bindings."""

from ._runtime_loader import preload_matching_netcdf_c

preload_matching_netcdf_c()

from . import clubb_api
from . import derived_types
from . import CLUBB_core

__all__ = ["clubb_api", "derived_types", "CLUBB_core"]
