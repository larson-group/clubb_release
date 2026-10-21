"""Converters for ConfigFlags <-> Fortran module storage."""

from typing import Iterable, Tuple

import clubb_f2py

from clubb_python.derived_types.config_flags import ConfigFlags, CONFIG_FLAG_NAMES

_NUM_INT_FLAGS = 9


def decode_config_flags(values: Iterable[object]) -> ConfigFlags:
    """Decode Fortran config flag sequence into ConfigFlags."""
    vals = list(values)
    if len(vals) != len(CONFIG_FLAG_NAMES):
        raise ValueError(
            f"Expected {len(CONFIG_FLAG_NAMES)} config flags, got {len(vals)}"
        )

    decoded = []
    for i, value in enumerate(vals):
        if i < _NUM_INT_FLAGS:
            decoded.append(int(value))
        else:
            decoded.append(bool(int(value)))
    return ConfigFlags(*decoded)


def encode_config_flags(flags: ConfigFlags) -> Tuple[object, ...]:
    """Encode ConfigFlags into Fortran args (ints + logicals)."""
    encoded = []
    for i, value in enumerate(flags):
        if i < _NUM_INT_FLAGS:
            encoded.append(int(value))
        else:
            encoded.append(bool(value))
    return tuple(encoded)


def set_fortran_config_flags(flags: ConfigFlags):
    """Push ConfigFlags into Fortran module storage."""
    clubb_f2py.f2py_init_config_flags(*encode_config_flags(flags))
