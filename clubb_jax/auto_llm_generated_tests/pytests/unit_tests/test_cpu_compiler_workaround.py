"""Check CPU compiler policy and the array-slice regression from JAX #42185.

Run through human_reviewed_tests/run_pytests.sh -jax. The numerical check
uses a fresh CPU process so XLA reads the model-import policy before compiling;
it requires no SCM case, benchmark files or Fortran library.
"""

import os
import shlex
import subprocess
import sys
import textwrap

import pytest

from clubb_jax.backends import cpu
from utilities.output_paths import REPO_ROOT


@pytest.mark.parametrize("flags", [
    "",
    '--xla_dump_to="/tmp/clubb dump" --xla_cpu_enable_fast_math=false',
    "--xla_backend_extra_options=--slp-vectorize-hor=false",
    "--xla_backend_extra_options=--combiner-disabled=false,--slp-vectorize-hor=false",
    "--xla_backend_extra_options '--combiner-disabled=false,--slp-vectorize-hor=false'",
])
def test_cpu_policy_preserves_other_flags_and_enforces_workaround(flags):
    """Keep inherited XLA/LLVM options while enforcing an idempotent CPU remedy."""
    env = {"XLA_FLAGS": flags}
    cpu.configure_environment(env)
    configured = env["XLA_FLAGS"]
    tokens = shlex.split(configured)
    extra, = [token.split("=", 1)[1] for token in tokens
              if token.startswith("--xla_backend_extra_options=")]
    options = dict(option.split("=", 1) for option in extra.split(","))
    assert options["--combiner-disabled"] == "true"
    if "slp-vectorize-hor" in flags:
        assert options["--slp-vectorize-hor"] == "false"
    if "xla_dump_to" in flags:
        assert "--xla_dump_to=/tmp/clubb dump" in tokens
        assert "--xla_cpu_enable_fast_math=false" in tokens
    cpu.configure_environment(env)
    assert env["XLA_FLAGS"] == configured


@pytest.mark.parametrize("flags", [
    "",
    "--xla_backend_extra_options=--combiner-disabled=false",
    '--xla_cpu_enable_fast_math=false --xla_dump_to="{dump_dir}"',
])
def test_direct_model_import_preserves_independent_rows_in_cpu_jit(flags, tmp_path):
    """Compare ordinary compiled array operations with NumPy after model import.

    This is the valid-index reproducer from https://github.com/jax-ml/jax/issues/42185.
    The affected compiler duplicates one row even without differentiation.
    """
    script = textwrap.dedent('''
        import clubb_jax.src
        import jax
        import jax.numpy as jnp
        import numpy as np

        jax.config.update("jax_enable_x64", True)

        def kernel(weight, read_index, clear, write_index, valid, a, b):
            x = a + b
            write_start = (write_index[0], jnp.int32(0), write_index[1])
            old_slice = jax.lax.dynamic_slice(x, write_start, (1, 2, 86))
            x = jax.lax.dynamic_update_slice(
                x, jnp.where(clear, 0.0, old_slice), write_start
            )
            selected = jax.lax.dynamic_slice(
                x, (read_index[0], jnp.int32(0), read_index[1]), (1, 2, 86)
            )
            result = jnp.where(valid, selected, 0.0).reshape(2, 86) * weight
            return jnp.pad(result, ((0, 0), (1, 1)))

        a = np.arange(352, dtype=np.float64).reshape(2, 2, 88)
        b = np.ones_like(a)
        weight = np.full((2, 86), 0.025)
        actual = np.asarray(jax.jit(kernel)(
            weight, np.array([0, 1], dtype=np.int32), np.array(True),
            np.array([1, 1], dtype=np.int32), np.array(True), a, b
        ))
        expected = np.pad((a + b)[0, :, 1:87] * weight, ((0, 0), (1, 1)))
        assert actual.shape == expected.shape == (2, 88)
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
    ''')
    env = dict(os.environ, CLUBB_JAX_ACCELERATOR="cpu", JAX_PLATFORMS="cpu",
               JAX_ENABLE_COMPILATION_CACHE="false")
    env["XLA_FLAGS"] = flags.format(dump_dir=tmp_path / "compiler dumps")
    result = subprocess.run([sys.executable, "-c", script], cwd=REPO_ROOT, env=env,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("accelerator", ["cuda13", "rocm", "metal"])
def test_gpu_model_import_leaves_compiler_flags_unchanged(accelerator):
    """Import model policy for a declared GPU without initializing that backend."""
    flags = "--xla_gpu_autotune_level=4"
    env = dict(os.environ, CLUBB_JAX_ACCELERATOR=accelerator, XLA_FLAGS=flags)
    result = subprocess.run(
        [sys.executable, "-c",
         'import os; import clubb_jax.src; print(os.environ["XLA_FLAGS"])'],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == flags
