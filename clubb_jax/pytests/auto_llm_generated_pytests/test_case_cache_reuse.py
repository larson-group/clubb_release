"""Check that fresh initialization retains JIT programs without stale physics.

Fortran's module-owned KK parameters remain authoritative. The JAX adaptation
must invalidate captured constants when inputs change, while output paths alone
must not invalidate an otherwise identical case.
"""

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.pytests.microphysics_test_inputs import initialize
from clubb_jax.src.Microphys import microphys_init_cleanup as lifecycle
from clubb_jax.src.Microphys.KK_microphys import parameters_KK
from clubb_jax.src.CLUBB_core.calc_pressure import init_pressure
from clubb_jax.src.CLUBB_core.grid_class import setup_grid


def test_repeated_configuration_reuses_trace_but_changed_constant_is_retraced(monkeypatch):
    """Reuse a captured KK constant for identical settings and retrace when it changes."""
    monkeypatch.setattr(lifecycle, '_initialization_settings', None)
    traces = []

    @jax.jit
    def read_constant():
        """Capture the module constant and record traces rather than executions."""
        traces.append(parameters_KK.C_evap)
        return jnp.asarray(parameters_KK.C_evap)

    try:
        initialize(microphys_scheme='khairoutdinov_kogan', c_evap=0.86, output_dir='first')
        assert float(read_constant()) == 0.86
        lifecycle.cleanup_microphys()
        initialize(microphys_scheme='khairoutdinov_kogan', c_evap=0.86, output_dir='second')
        assert float(read_constant()) == 0.86
        assert traces == [0.86]
        lifecycle.cleanup_microphys()
        initialize(microphys_scheme='khairoutdinov_kogan', c_evap=0.7, output_dir='third')
        assert float(read_constant()) == 0.7
        assert traces == [0.86, 0.7]
        lifecycle.cleanup_microphys()
        initialize(microphys_scheme='khairoutdinov_kogan')
        assert float(read_constant()) == 0.86
    finally:
        lifecycle.cleanup_microphys()


def test_pressure_jit_preserves_eager_initialization_values():
    """Compare two-column hydrostatic profiles from eager and compiled initialization."""
    gr = setup_grid(2, 40., 0., 400.)
    thvm = jnp.broadcast_to(jnp.linspace(290., 310., gr.nzt), (2, gr.nzt))
    p_sfc = jnp.asarray([100000., 99000.])
    expected = init_pressure.__wrapped__(thvm, p_sfc, gr)
    actual = init_pressure(thvm, p_sfc, gr)
    for value, reference in zip(actual, expected):
        assert value.shape == reference.shape
        np.testing.assert_allclose(value, reference, rtol=3e-15, atol=1e-10)
