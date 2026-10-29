"""Launcher cache policy and reuse across isolated processes/cache clearing."""

import json
import os
from pathlib import Path
import subprocess
import sys

from clubb_jax import run_jax


def test_launcher_cache_survives_checkout_location_and_respects_override(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.delenv("JAX_COMPILATION_CACHE_DIR", raising=False)
    expected = tmp_path / ".cache" / "clubb-jax" / "compilation" / "cpu"
    assert run_jax.runtime_configuration("cpu")["environment"]["JAX_COMPILATION_CACHE_DIR"] == str(expected)
    assert not expected.exists()  # Inspection does not create runtime files.
    monkeypatch.setenv("JAX_COMPILATION_CACHE_DIR", str(tmp_path / "custom"))
    assert run_jax.runtime_configuration("cpu")["environment"]["JAX_COMPILATION_CACHE_DIR"] == str(tmp_path / "custom")


def test_compiled_program_reuses_disk_cache_in_fresh_processes(tmp_path):
    code = """
import json
import jax
import jax.numpy as jnp

@jax.jit
def cached_sum(x):
    return jnp.sin(x).sum()

x = jnp.arange(4096, dtype=jnp.float32) * 0.01
first = float(cached_sum(x))
jax.clear_caches()
second = float(cached_sum(x))
print(json.dumps([first, second]))
"""
    env = os.environ | {
        "JAX_PLATFORMS": "cpu",
        "JAX_COMPILATION_CACHE_DIR": str(tmp_path / "cache"),
        "JAX_PERSISTENT_CACHE_MIN_COMPILE_TIME_SECS": "0",
        "JAX_LOG_COMPILES": "1",
    }
    cold = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
    warm = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
    assert json.loads(cold.stdout) == json.loads(warm.stdout)
    assert json.loads(warm.stdout)[0] == json.loads(warm.stdout)[1]
    assert any((tmp_path / "cache").iterdir())
    assert "Persistent compilation cache hit for 'jit_cached_sum'" in warm.stderr
