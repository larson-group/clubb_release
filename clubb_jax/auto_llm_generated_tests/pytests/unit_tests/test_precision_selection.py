"""CLUBB_JAX_PRECISION parsing, launcher reporting and the package-wide JAX setting."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from clubb_jax import run_jax
from clubb_jax.backends.common import precision_name

REPO_ROOT = Path(__file__).resolve().parents[4]


@pytest.mark.parametrize("value, expected", [
    (None, "double"), ("", "double"), ("double", "double"), ("float64", "double"),
    ("single", "single"), (" SP ", "single"), ("float32", "single"), ("real4", "single"),
])
def test_precision_name_aliases(value, expected):
    assert precision_name(value) == expected


def test_launcher_reports_selected_precision(monkeypatch):
    monkeypatch.delenv("CLUBB_JAX_PRECISION", raising=False)
    assert run_jax.runtime_configuration("cpu")["precision"] == "double"
    monkeypatch.setenv("CLUBB_JAX_PRECISION", "f32")
    assert run_jax.runtime_configuration("cpu")["precision"] == "single"


@pytest.mark.parametrize("value, x64, dtype", [("double", True, "float64"), ("single", False, "float32")])
def test_importing_model_package_applies_precision(value, x64, dtype):
    # A fresh process imports one model module directly; no module calls the setup itself.
    code = """
import json
import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.new_pdf import calc_coef_wp4_implicit
coef = calc_coef_wp4_implicit(jnp.full(3, 0.3), jnp.full(3, 0.5), jnp.full(3, 0.4), jnp.full(3, 0.6))
print(json.dumps({"x64": jax.config.jax_enable_x64, "dtype": str(coef.dtype)}))
"""
    env = {**os.environ, "CLUBB_JAX_PRECISION": value, "JAX_PLATFORMS": "cpu"}
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=REPO_ROOT, env=env,
        capture_output=True, text=True, check=True,
    )
    assert json.loads(result.stdout.strip().splitlines()[-1]) == {"x64": x64, "dtype": dtype}
