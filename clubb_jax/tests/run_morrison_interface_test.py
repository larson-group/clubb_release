#!/usr/bin/env python3
"""Validate eager/JIT Morrison species interfaces and compiled fatal diagnostics."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.parse_args()
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()

from contextlib import redirect_stdout
import io
from types import SimpleNamespace
from unittest.mock import patch
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
from clubb_jax.src.Microphys.morrison_microphys_module import morrison_microphys_driver
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.Microphys.microphys_init_cleanup import cleanup_microphys
from clubb_jax.pytests.microphysics_test_inputs import initialize


def check_morrison_species_interface_eager_and_jit(dimension, ncol):
    _, _, metadata, *_ = initialize(microphys_scheme='morrison',
        l_ice_microphys=dimension>2, l_graupel=dimension>6)
    nzt=5
    gr = SimpleNamespace(zt=jnp.broadcast_to(jnp.arange(nzt)*100.,(ncol,nzt)),
                         dzt=jnp.full((ncol,nzt),100.))
    one=jnp.ones((ncol,nzt)); zero=jnp.zeros_like(one)
    hydromet=jnp.zeros((ncol,nzt,dimension))
    for i in range(0,dimension,2):
        hydromet=hydromet.at[...,i].set(1.e-5)
        hydromet=hydromet.at[...,i+1].set(2.e4)
    # Columns differ to catch accidental flattening or cross-column transport.
    hydromet=hydromet*jnp.arange(1,ncol+1)[:,None,None]
    stats=JaxStats.empty(l_sample=True, names=("rrm_auto", "precip_rate_sfc"),
        grids=("zt", "sfc"), ncol=ncol, max_nlev=nzt)
    def run(hm):
        return morrison_microphys_driver(
            gr, ncol, 10., nzt,                                                # In
            dimension, metadata,                                               # In
            False, jnp.linspace(258.,285.,nzt)[None,:]*one, zero, 90000.*one,  # In
            one, one, 0.5*one, 0.2*one,                                        # In
            100.*one, 1.e-4*one, 1.e8*one, zero, 0.007*one, hm,                # In
            1,                                                                 # In
            one,                                                               # In
            stats,                                                             # InOut
        )
    eager=run(hydromet)
    compiled=jax.jit(run)(hydromet)
    assert len(eager)==12
    for a,b in zip(eager[1:],compiled[1:]):
        assert np.all(np.isfinite(a))
        np.testing.assert_allclose(a,b,rtol=1.e-6,atol=1.e-9)
    assert eager[1].shape==hydromet.shape
    assert np.all(np.asarray(hydromet+10.*eager[1])>=-1.e-12)
    assert np.all(np.asarray(eager[2][...,1:])==0.)
    assert np.all(np.asarray(eager[2][...,0])<=0.)
    for a,b in zip(eager[0].buffers,compiled[0].buffers):
        np.testing.assert_allclose(a,b,rtol=1.e-6,atol=1.e-9)
    np.testing.assert_allclose(eager[0].buffers[0][0],eager[7])
    assert np.all(np.asarray(eager[0].nsamples[0]) == 1)
    assert np.all(np.asarray(eager[0].nsamples[2]) == 1)


def check_morrison_debug_output_is_live_in_jit():
    from clubb_jax.src.CLUBB_core import error_code
    from clubb_jax.src.Microphys import morrison_microphys_module as morrison
    _, _, metadata, *_ = initialize(microphys_scheme='morrison')
    core = morrison.M2005MICRO_GRAUPEL
    def broken_core(*args):
        output = core(*args)
        output['T3DTEN'] = jnp.full_like(output['T3DTEN'], jnp.nan)
        return output
    one = jnp.ones((1, 3)); zero = jnp.zeros_like(one)
    gr = SimpleNamespace(zt=100*one, dzt=100*one)
    stats = JaxStats.empty(l_sample=False, names=(), ncol=1, max_nlev=3)
    def run(hydromet):
        return morrison.morrison_microphys_driver(
            gr, 1, 10., 3,                                           # In
            2, metadata,                                             # In
            False, 280*one, zero, 90000*one,                         # In
            one, one, .5*one, .2*one,                                # In
            100*one, 1.e-4*one, 1.e8*one, zero, .007*one, hydromet,  # In
            1,                                                       # In
            one,                                                     # In
            stats,                                                   # InOut
        )
    output = io.StringIO()
    # Wait for host callbacks before inspecting and restoring the injected diagnostic.
    with patch.object(error_code, '_debug_level', 2), \
            patch.object(morrison, 'M2005MICRO_GRAUPEL', broken_core), \
            redirect_stdout(output):
        jax.block_until_ready(jax.jit(run)(jnp.zeros((1, 3, 2))))
        jax.effects_barrier()
    assert 'non-finite detected in a Morrison microphysics tendency' in output.getvalue()


if __name__ == "__main__":
    for dimension, ncol in ((2, 1), (6, 2), (8, 2)):
        try:
            check_morrison_species_interface_eager_and_jit(dimension, ncol)
            print(f"Morrison interface passed: {dimension} species, {ncol} columns.")
        finally:
            cleanup_microphys()
            jax.clear_caches()
    try:
        check_morrison_debug_output_is_live_in_jit()
        print("Morrison compiled diagnostic passed.")
    finally:
        cleanup_microphys()
        jax.clear_caches()
