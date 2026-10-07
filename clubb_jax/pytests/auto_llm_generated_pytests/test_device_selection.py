"""Device requests cross the public launcher boundary before JAX initializes.

Exercise inspection and execution with a synthetic CUDA inventory, including
explicit false overriding an inherited true setting. No GPU or JAX import is
needed for the launcher to validate and prepare each child's environment.
"""

import os

import pytest

from clubb_jax import backends, run_jax


GPU_A = "GPU-aaaaaaaa-1111-2222-3333-000000000001"
GPU_B = "GPU-bbbbbbbb-1111-2222-3333-000000000002"


@pytest.mark.parametrize("prealloc_gpu_mem", [None, False, True])
def test_selection_reaches_inspection_and_execution_without_mutating_parent(monkeypatch, tmp_path, prealloc_gpu_mem):
    monkeypatch.setattr(backends, "native_gpu_accelerator", lambda: "cuda13")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", GPU_A)
    monkeypatch.setenv("XLA_PYTHON_CLIENT_PREALLOCATE", "true")
    monkeypatch.setattr(backends.cuda, "query_gpus", lambda: ([
        {"uuid": GPU_A, "index": "0", "pci_bus_id": "0000:01:00.0"},
        {"uuid": GPU_B, "index": "1", "pci_bus_id": "0000:02:00.0"},
    ], ""))
    inspected = []
    executed = []

    def inspect(*args, **kwargs):
        inspected.append(args[-1].copy())
        return 0

    monkeypatch.setattr(run_jax, "_run_inspection", inspect)
    monkeypatch.setattr(run_jax, "_prepare_environment", lambda *args: tmp_path / "python")
    monkeypatch.setattr(run_jax, "_print_runtime_summary", lambda *args: None)
    monkeypatch.setattr(run_jax.os, "chdir", lambda *args: None)
    monkeypatch.setattr(run_jax.os, "execvpe", lambda python, args, env: executed.append((args, env.copy())))
    args = run_jax.runtime_arguments("gpu", device=GPU_B, prealloc_gpu_mem=prealloc_gpu_mem)
    assert run_jax.main([*args, "-info=json"]) == 0
    assert run_jax.main([*args, "case.in"]) == 0
    expected_preallocation = "true" if prealloc_gpu_mem is None else str(prealloc_gpu_mem).lower()
    for env in [*inspected, executed[0][1]]:
        assert env["CUDA_VISIBLE_DEVICES"] == GPU_B
        assert env["XLA_PYTHON_CLIENT_PREALLOCATE"] == expected_preallocation
        assert env["CLUBB_JAX_PROFILE"] == "gpu"
    assert executed[0][0][-1] == "case.in"
    assert os.environ["CUDA_VISIBLE_DEVICES"] == GPU_A
    assert os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] == "true"


@pytest.mark.parametrize("device", ["1", "GPU-short", GPU_A + "," + GPU_B, "MIG-example", "bad\nvalue"])
def test_public_selection_rejects_invalid_device_identifiers(device):
    with pytest.raises(ValueError):
        run_jax.runtime_selection("gpu", device=device)


@pytest.mark.parametrize("args", [
    ["-options=cpu,device=" + GPU_A],
    ["-options=gpu,xla_prealloc,prealloc_gpu_mem=false"],
    ["-options=gpu,prealloc_gpu_mem=maybe"],
    [f"-options=gpu,device={GPU_A},device={GPU_B}"],
    ["-options=rocm,device=" + GPU_A],
])
def test_invalid_launch_selection_fails_before_setup(monkeypatch, args):
    monkeypatch.setattr(run_jax, "_prepare_environment", lambda *args: pytest.fail("invalid selection reached setup"))
    with pytest.raises((ValueError, run_jax.LauncherError)):
        run_jax.main(args)


def test_default_selection_preserves_visibility_and_configuration_is_read_only(monkeypatch):
    monkeypatch.setattr(backends, "native_gpu_accelerator", lambda: "cuda13")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", GPU_A)
    monkeypatch.setenv("XLA_PYTHON_CLIENT_PREALLOCATE", "true")
    assert run_jax.runtime_configuration("gpu")["environment"]["CUDA_VISIBLE_DEVICES"] == GPU_A
    configuration = run_jax.runtime_configuration("gpu", device=GPU_B, prealloc_gpu_mem=False)
    assert configuration["environment"]["CUDA_VISIBLE_DEVICES"] == GPU_B
    assert configuration["environment"]["XLA_PYTHON_CLIENT_PREALLOCATE"] == "false"
    assert os.environ["CUDA_VISIBLE_DEVICES"] == GPU_A
    assert os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] == "true"


def test_scm_wrappers_forward_the_selection_to_the_launcher():
    from types import SimpleNamespace
    from run_scripts.run_scm import extract_jax_options
    from run_scripts.run_scm import choose_run_command

    args = run_jax.runtime_arguments("gpu", device=GPU_B, prealloc_gpu_mem=False, scm=True)
    assert args == [f"-jax=gpu,device={GPU_B},prealloc_gpu_mem=false"]
    normalized, value, occurrences = extract_jax_options([*args, "bomex"])
    assert normalized == ["-jax", "bomex"] and occurrences == 1
    command, _, _ = choose_run_command(SimpleNamespace(
        exe=None, python=False, jax=True, jax_options=value, gdb=False,
    ))
    assert command[1:] == ["-options=" + value]
    values, remaining = run_jax.parse_launcher_args(["-options=" + value, "case.in"])
    assert values["device"] == GPU_B
    assert values["prealloc_gpu_mem"] is False
    assert remaining == ["case.in"]


@pytest.mark.parametrize("memory", [False, True])
def test_encoded_settings_round_trip_through_the_public_selection_interface(memory):
    value = f"gpu,device={GPU_B},prealloc_gpu_mem={str(memory).lower()}"
    selection = run_jax.runtime_selection(value)
    assert selection == {"profile": "gpu", "device": GPU_B, "prealloc_gpu_mem": memory}
    assert run_jax.runtime_arguments(value, scm=True) == ["-jax=" + value]
    assert run_jax.runtime_selection(value, device=GPU_B, prealloc_gpu_mem=memory) == selection
    with pytest.raises(ValueError, match="conflicts"):
        run_jax.runtime_selection(value, device=GPU_A)
    with pytest.raises(ValueError, match="conflicts"):
        run_jax.runtime_selection(value, prealloc_gpu_mem=not memory)


def test_legacy_memory_modifier_remains_compatible():
    assert run_jax.runtime_selection("gpu,xla_prealloc") == run_jax.runtime_selection("gpu,prealloc_gpu_mem=true")


@pytest.mark.parametrize("value", [
    "gpu,device", "gpu,prealloc_gpu_mem", "gpu,preallocation=false",
    "gpu,prealloc_gpu_mem=false,prealloc_gpu_mem=true",
])
def test_malformed_selection_modifiers_fail_in_the_launcher(value):
    with pytest.raises(ValueError):
        run_jax.runtime_selection(value)
