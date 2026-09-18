"""Per-job JAX device selection, independent of CUDA's numeric device order."""

import os
import re
import shlex


GPU_UUID_PATTERN = r"^(?:GPU-[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12})?$"


def normalize_jax_gpu(value):
    """An empty selection inherits the launcher environment; otherwise use a UUID."""
    value = str(value or "").strip()
    if not re.fullmatch(GPU_UUID_PATTERN, value):
        raise ValueError("JAX GPU must be a full GPU UUID or empty for Default")
    return value


def jax_device_env(settings, base_env=None):
    """Copy the environment and apply only this job's explicit GPU selection."""
    env = dict(os.environ if base_env is None else base_env)
    prealloc = settings.get("jax_xla_prealloc")
    if prealloc is not None:
        if not isinstance(prealloc, bool):
            raise ValueError("jax_xla_prealloc must be a boolean")
        if settings.get("implementation") != "jax" or settings.get("jax_profile") != "gpu":
            if prealloc:
                raise ValueError("jax_xla_prealloc requires the JAX GPU profile")
        else:
            env["XLA_PYTHON_CLIENT_PREALLOCATE"] = str(prealloc).lower()
    gpu = normalize_jax_gpu(settings.get("jax_gpu"))
    if gpu:
        if settings.get("implementation") != "jax" or settings.get("jax_profile") != "gpu":
            raise ValueError("jax_gpu requires the JAX GPU profile")
        env["CUDA_VISIBLE_DEVICES"] = gpu
    return env


def jax_command_display(command, settings):
    """Include explicit per-job environment settings in copyable commands."""
    overrides = jax_device_env(settings, {})
    prefix = " ".join(f"{key}={shlex.quote(value)}" for key, value in overrides.items())
    return (prefix + " " if prefix else "") + shlex.join(command)
