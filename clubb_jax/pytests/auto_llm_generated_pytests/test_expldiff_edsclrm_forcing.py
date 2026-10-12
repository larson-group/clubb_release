"""With l_do_expldiff_rtm_thlm, the thlm/rtm eddy-scalar copies get no forcing."""

from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.time_dependent_input import (
    apply_time_dependent_forcings_from_array,
)
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
from clubb_jax.src.Input_fields.input_reader import two_dim_read_var
from clubb_jax.src.Input_fields.input_names import rt_f_name, thetal_f_name


@pytest.mark.parametrize("l_do_expldiff_rtm_thlm", [False, True])
def test_expldiff_copies_get_no_time_dependent_forcing(monkeypatch, l_do_expldiff_rtm_thlm):
    nzt, ngrdcol, edsclr_dim = 4, 2, 3
    names = [thetal_f_name, rt_f_name]
    names += [rt_f_name] * (time_dependent_input.nforcings - 1 - len(names))
    forcing_data = [two_dim_read_var("z[m]", "Time[s]", "z[m]", np.zeros((nzt, 1)))]
    forcing_data += [two_dim_read_var(name, "Time[s]", "z[m]", np.zeros((nzt, 1)))
                     for name in names]
    monkeypatch.setattr(time_dependent_input, "t_dependent_forcing_data", forcing_data)

    thl_forcing = jnp.linspace(-2.0e-5, 1.0e-5, nzt)
    rt_forcing = jnp.linspace(-4.0e-8, 1.0e-8, nzt)
    # Rows after thetal_f and rtm_f are -999.9, i.e. missing, and are skipped.
    forcings_array = jnp.full((time_dependent_input.nforcings, nzt), -999.9)
    forcings_array = forcings_array.at[1].set(thl_forcing).at[2].set(rt_forcing)

    sclr_idx = SclrIdx(iisclr_rt=-1, iisclr_thl=-1, iisclr_CO2=-1,
                       iiedsclr_rt=3, iiedsclr_thl=2, iiedsclr_CO2=1)
    zeros_zt = jnp.zeros((ngrdcol, nzt))
    edsclrm_forcing = jnp.full((ngrdcol, nzt, edsclr_dim), 7.0)

    result = apply_time_dependent_forcings_from_array(
        ngrdcol, nzt - 1, nzt, 0, edsclr_dim, sclr_idx, l_do_expldiff_rtm_thlm,
        SimpleNamespace(nzt=nzt, nzm=nzt - 1, ngrdcol=ngrdcol),
        zeros_zt, jnp.ones_like(zeros_zt), jnp.ones_like(zeros_zt),
        forcings_array,
        zeros_zt, zeros_zt, zeros_zt, zeros_zt, zeros_zt, zeros_zt,
        zeros_zt, jnp.zeros((ngrdcol, nzt - 1)), zeros_zt, zeros_zt,
        jnp.zeros((ngrdcol, nzt, 0)), edsclrm_forcing,
    )
    thlm_f, rtm_f, edsclrm_forcing = result[0], result[1], np.asarray(result[-1])

    # thlm and rtm themselves are always forced.
    np.testing.assert_array_equal(thlm_f, jnp.broadcast_to(thl_forcing, (ngrdcol, nzt)))
    np.testing.assert_array_equal(rtm_f, jnp.broadcast_to(rt_forcing, (ngrdcol, nzt)))
    # Other eddy scalars keep their forcing.
    np.testing.assert_array_equal(edsclrm_forcing[:, :, 0], 7.0)
    if l_do_expldiff_rtm_thlm:
        np.testing.assert_array_equal(edsclrm_forcing[:, :, 1], 0.0)
        np.testing.assert_array_equal(edsclrm_forcing[:, :, 2], 0.0)
    else:
        np.testing.assert_array_equal(edsclrm_forcing[:, :, 1], thlm_f)
        np.testing.assert_array_equal(edsclrm_forcing[:, :, 2], rtm_f)
