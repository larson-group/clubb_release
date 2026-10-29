"""Invalid NetCDF fields preserve the restart reader's source error contract."""

import jax.numpy as jnp
from netCDF4 import Dataset
import numpy as np
import pytest

from clubb_jax.src.Input_fields import input_fields


@pytest.mark.parametrize('kind,units', [('i4', 'K'), ('i8', 'K'), ('f8', None), ('f8', 1)])
def test_invalid_field_returns_read_error_without_replacing_state(tmp_path, kind, units):
    path = tmp_path / 'restart.nc'
    with Dataset(path, 'w') as ds:
        for dim, size in [('time', 1), ('zt', 2), ('col', 1)]:
            ds.createDimension(dim, size)
        ds.createVariable('zt', 'f8', ('zt',))[:] = [0., 10.]
        var = ds.createVariable('thlm', kind, ('time', 'zt', 'col'))
        var[:] = np.full((1, 2, 1), 300)
        if units is not None:
            var.units = units
    original = jnp.array([[280., 281.]])
    result, error = input_fields.get_clubb_variable_interpolated(
        True, path, 'thlm', 2, 1, np.array([0., 10.]), original,
    )
    assert error
    np.testing.assert_array_equal(result, original)
