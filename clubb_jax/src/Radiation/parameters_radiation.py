"""JAX-owned equivalent of ``Radiation/parameters_radiation.F90``.

Description:
    Parameters for radiation schemes

References:
    None

JAX adaptation:
    The source stores these fields in mutable module state. A compiled JAX
    calculation instead receives immutable ``RadiationParameters`` metadata.
    Source lookup tables are dynamic pytree leaves; scheme selection and scalar
    namelist values are static. ``l_soil_veg`` and LBA tables are explicit
    additions from their respective source modules. BUGSrad-only source fields,
    including ``sol_const``, are excluded because that configuration is rejected
    at initialization.

Excluded source field comments:
    ``sol_const`` - Solar constant
    ``radiation_top`` - The top of the atmosphere fed into a radiation scheme.
                   The computational grid should be extended to reach this
                   altitude.
    ``alndr`` - Near-IR direct surface albedo   [-]
    ``alvdf`` - Visible diffuse surface albedo  [-]
    ``alndf`` - Near-IR diffuse surface albedo  [-]
    ``slr`` - Fraction of daylight
    ``l_use_default_std_atmosphere``:
        Flag to signal the use of the U.S. Standard Atmosphere Profile, 1976
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np


_RADIATION_LIST_SIZE = 20


# JAX host setup must normalize source namelist lists before they become fixed
# shape compiled parameters; Fortran assigns those module arrays directly.
def _float_array(value, length: int, name: str):
    """Normalize a fixed-size source radiation list at the host boundary."""
    array = np.asarray(value if isinstance(value, (list, tuple, np.ndarray)) else [value], dtype=np.float64)
    if array.size > length:
        raise ValueError(f"{name} has {array.size} values; parameters_radiation allows {length}.")
    return jnp.asarray(np.pad(array, (0, length - array.size)))


@jax.tree_util.register_pytree_node_class
@dataclass(frozen=True)
class RadiationParameters:
    """Immutable replacement for active ``parameters_radiation`` fields.

    Retained fields through ``nparam`` follow the corresponding source module
    declaration order. The final fields are documented functional inputs from
    the soil and LBA source modules.
    """

    rad_scheme: str  # Either BUGSrad, simplified, or simplied_bomex

    # Albedo values (alvdr is used in the simplifed schemes as well)
    alvdr: float  # Visible direct surface albedo   [-]

    # Long-wave constants (simplified radiation)
    kappa: float  # A constant (Duynkerke eqn. 5)                   [m^2/kg]
    F0: float  # Coefficient for cloud top heating (see Stevens) [W/m^2]
    F1: float  # Coefficient for cloud base heating (see Stevens)[W/m^2]

    # Short-wave constants
    eff_drop_radius: float  # Effective droplet radius [m]
    gc: float  # Asymmetry parameter, "g" in Duynkerke           [-]
    omega: float  # Single-scattering albedo                        [-]

    Fs_values: object  # List of Fs0 values for simplified radiation
    cos_solar_zen_times: object  # List of cosine of the solar zenith angle times
    cos_solar_zen_values: object  # List of cosine of the solar zenith angle values
    l_fix_cos_solar_zen: bool
    l_sw_radiation: bool
    l_rad_above_cloud: bool  # Use DYCOMS II RF02 heaviside step function
    nparam: int
    l_soil_veg: bool
    lba_zrad: object  # Altitudes        [m]
    lba_krad: object  # Radiative tendencies     [K/s]

    def tree_flatten(self):
        children = (
            self.Fs_values,
            self.cos_solar_zen_times,
            self.cos_solar_zen_values,
            self.lba_zrad,
            self.lba_krad,
        )
        aux_data = (
            self.rad_scheme, self.alvdr, self.kappa, self.F0, self.F1,
            self.eff_drop_radius, self.gc, self.omega,
            self.l_fix_cos_solar_zen, self.l_sw_radiation,
            self.l_rad_above_cloud, self.nparam, self.l_soil_veg,
        )
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        (
            rad_scheme, alvdr, kappa, F0, F1, eff_drop_radius, gc, omega,
            l_fix_cos_solar_zen, l_sw_radiation, l_rad_above_cloud, nparam,
            l_soil_veg,
        ) = aux_data
        (
            Fs_values, cos_solar_zen_times, cos_solar_zen_values,
            lba_zrad, lba_krad,
        ) = children
        return cls(
            rad_scheme=rad_scheme,
            alvdr=alvdr,
            kappa=kappa,
            F0=F0,
            F1=F1,
            eff_drop_radius=eff_drop_radius,
            gc=gc,
            omega=omega,
            Fs_values=Fs_values,
            cos_solar_zen_times=cos_solar_zen_times,
            cos_solar_zen_values=cos_solar_zen_values,
            l_fix_cos_solar_zen=l_fix_cos_solar_zen,
            l_sw_radiation=l_sw_radiation,
            l_rad_above_cloud=l_rad_above_cloud,
            nparam=nparam,
            l_soil_veg=l_soil_veg,
            lba_zrad=lba_zrad,
            lba_krad=lba_krad,
        )


def initialize_radiation_parameters(cfg: dict, rad_scheme: str, case_setups_dir: Path) -> RadiationParameters:
    """Construct immutable radiation-module values at the host I/O boundary.

    This target-only constructor replaces source namelist/module assignment.
    ``cfg`` supplies source parameter values in declaration order; LBA data is
    read only when the source ``rad_scheme`` selects LBA.
    """
    from clubb_jax.src.Radiation.simple_rad_module import simple_rad_lba_init

    case_setups_dir = Path(case_setups_dir)
    Fs_values = _float_array(cfg.get("fs_values", 0.0), _RADIATION_LIST_SIZE, "Fs_values")
    cos_solar_zen_times = _float_array(
        cfg.get("cos_solar_zen_times", 0.0), _RADIATION_LIST_SIZE, "cos_solar_zen_times",
    )
    cos_solar_zen_values = _float_array(
        cfg.get("cos_solar_zen_values", 0.0), _RADIATION_LIST_SIZE, "cos_solar_zen_values",
    )
    nparam = int(np.asarray(cfg.get("cos_solar_zen_values", [0.0])).size)

    if bool(cfg.get("l_sw_radiation", False)) and not bool(cfg.get("l_fix_cos_solar_zen", False)):
        day = int(cfg.get("day", 0))
        month = int(cfg.get("month", 0))
        year = int(cfg.get("year", 0))
        if not 1 <= month <= 12:
            raise ValueError("Problem with Julian day conversion in gregorian2julian_day.")
        julian_day = day + sum((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[:month - 1])
        julian_day += int(year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) and month > 2)
        days_in_year = 366 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 365
        # Source ``gregorian2julian_day`` uses ``error stop`` for this condition.
        if julian_day > days_in_year:
            raise ValueError("Problem with Julian day conversion in gregorian2julian_day.")

    if rad_scheme == "lba":
        lba_zrad, lba_krad = simple_rad_lba_init(case_setups_dir / "lba_forcings")
    else:
        lba_zrad = jnp.zeros((33,), dtype=jnp.float64)
        lba_krad = jnp.zeros((33, 36), dtype=jnp.float64)

    return RadiationParameters(
        rad_scheme=rad_scheme,
        alvdr=float(cfg.get("alvdr", 0.0)),
        kappa=float(cfg.get("kappa", 0.0)),
        F0=float(cfg.get("f0", 0.0)),
        F1=float(cfg.get("f1", 0.0)),
        eff_drop_radius=float(cfg.get("eff_drop_radius", 10.0e-6)),
        gc=float(cfg.get("gc", 0.85)),
        omega=float(cfg.get("omega", 0.999)),
        Fs_values=Fs_values,
        cos_solar_zen_times=cos_solar_zen_times,
        cos_solar_zen_values=cos_solar_zen_values,
        l_fix_cos_solar_zen=bool(cfg.get("l_fix_cos_solar_zen", False)),
        l_sw_radiation=bool(cfg.get("l_sw_radiation", False)),
        l_rad_above_cloud=bool(cfg.get("l_rad_above_cloud", False)),
        nparam=nparam,
        l_soil_veg=bool(cfg.get("l_soil_veg", False)),
        lba_zrad=jnp.asarray(lba_zrad),
        lba_krad=jnp.asarray(lba_krad),
    )


__all__ = ["RadiationParameters", "initialize_radiation_parameters"]
