"""Tests for JAX runtime inspection and launcher preflight."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from clubb_jax import runtime_info


def test_cuda13_rejects_an_old_driver():
    selectable, reason = runtime_info._gpu_compatibility(
        [
            {
                "name": "Test GPU",
                "driver_version": "470.239",
                "compute_capability": "8.0",
            }
        ],
        "",
    )

    assert selectable is False
    assert "driver 580 or newer" in reason
    assert "470.239" in reason


def test_cuda13_accepts_a_three_component_new_driver_version():
    selectable, reason = runtime_info._gpu_compatibility(
        [
            {
                "name": "Test GPU",
                "driver_version": "610.57.04",
                "compute_capability": "8.0",
            }
        ],
        "",
    )

    assert selectable is True
    assert reason == ""


def test_cuda13_rejects_an_old_compute_capability():
    selectable, reason = runtime_info._gpu_compatibility(
        [
            {
                "name": "Test GPU",
                "driver_version": "580.10",
                "compute_capability": "7.0",
            }
        ],
        "",
    )

    assert selectable is False
    assert "compute capability 7.5 or newer" in reason


def test_gpu_query_retries_without_compute_capability(monkeypatch):
    responses = iter(
        (
            subprocess.CompletedProcess(
                [], 2, "", 'Field "compute_cap" is not a valid field to query.\n'
            ),
            subprocess.CompletedProcess(
                [], 0, "0, GPU-test, Test GPU, 580.10, 24576, 00000000:01:00.0\n", ""
            ),
        )
    )
    monkeypatch.setattr(runtime_info.shutil, "which", lambda _name: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(runtime_info.subprocess, "run", lambda *_args, **_kwargs: next(responses))

    gpus, error = runtime_info._query_nvidia_gpus()

    assert error == ""
    assert gpus == [
        {
            "index": "0",
            "uuid": "GPU-test",
            "name": "Test GPU",
            "driver_version": "580.10",
            "memory_mib": 24576,
            "compute_capability": "",
            "pci_bus_id": "00000000:01:00.0",
        }
    ]


def test_gpu_report_labels_expected_device_without_claiming_backend_ready(monkeypatch, tmp_path):
    monkeypatch.setenv("CUDA_DEVICE_ORDER", "PCI_BUS_ID")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1")
    gpus = [
        {
            "index": "0",
            "uuid": "GPU-test-a",
            "name": "Test GPU A",
            "driver_version": "610.57.04",
            "memory_mib": 6144,
            "compute_capability": "7.5",
            "pci_bus_id": "00000000:01:00.0",
        },
        {
            "index": "1",
            "uuid": "GPU-test-b",
            "name": "Test GPU B",
            "driver_version": "610.57.04",
            "memory_mib": 10240,
            "compute_capability": "8.6",
            "pci_bus_id": "00000000:02:00.0",
        },
    ]
    monkeypatch.setattr(runtime_info, '_query_nvidia_gpus', lambda: (gpus, ''))
    requirements = tmp_path / 'requirements.txt'
    requirements.touch()
    info = runtime_info.inspect_runtime(
        'gpu', 'cuda13', requirements, tmp_path / 'venv', '0.11.0'
    )

    output = runtime_info.format_human(info)

    assert " CLUBB JAX GPU PRE-RUN CHECK" in output
    assert " DEVICE SELECTION" in output
    assert " PHYSICAL GPU INVENTORY" in output
    assert "   0    --  hidden        Test GPU A" in output
    assert "   1     0  SELECTED      Test GPU B" in output
    assert " PLANNED DEVICE" in output
    assert " Expected JAX GPU 0 --> physical GPU 1 (Test GPU B)" in output
    assert "Backend initialization: not checked" in output
    assert info['runtime']['backend_verified'] is False


def test_numeric_visible_device_resolves_to_its_nvidia_smi_uuid(monkeypatch):
    monkeypatch.setattr(
        runtime_info,
        "_query_nvidia_gpus",
        lambda: (
            [
                {"index": "0", "uuid": "GPU-test-a"},
                {"index": "1", "uuid": "GPU-test-b"},
            ],
            "",
        ),
    )

    assert runtime_info.resolve_visible_devices("1,0") == "GPU-test-b,GPU-test-a"
    assert runtime_info.resolve_visible_devices("GPU-test-b") == "GPU-test-b"


def test_missing_environment_requires_setup_but_remains_selectable(tmp_path, monkeypatch):
    requirements = tmp_path / "requirements.txt"
    requirements.write_text("jax[cpu]==0.11.0\n", encoding="utf-8")
    monkeypatch.setattr(
        runtime_info,
        "_cpu_info",
        lambda: {"model": "Test CPU", "logical_cpus": 8},
    )

    info = runtime_info.inspect_runtime(
        "cpu",
        "cpu",
        requirements,
        tmp_path / "missing-venv",
        "0.11.0",
        "3.12",
    )

    assert info["selectable"] is True
    assert info["status"] == "setup_required"
    assert info["runtime"]["jax"]["required"] == "0.11.0"
    assert info["runtime"]["python"]["installed"] is False
    assert info["runtime"]["python"]["version"] == "3.12"


def test_human_output_labels_planned_versions(tmp_path, monkeypatch):
    requirements = tmp_path / "requirements.txt"
    requirements.write_text("jax[cpu]==0.11.0\n", encoding="utf-8")
    monkeypatch.setattr(
        runtime_info,
        "_cpu_info",
        lambda: {"model": "Test CPU", "logical_cpus": 8},
    )
    info = runtime_info.inspect_runtime(
        "cpu", "cpu", requirements, Path(tmp_path / "venv"), "0.11.0", "3.12"
    )

    output = runtime_info.format_human(info)

    assert " CLUBB JAX CPU PRE-RUN CHECK" in output
    assert " COMPUTE DEVICE" in output
    assert " CPU: Test CPU (8 logical CPUs)" in output
    assert " JAX: 0.11.0 (planned)" in output


def test_gpu_preflight_fails_before_creating_an_environment(tmp_path):
    repo_root = Path(__file__).resolve().parents[2]
    wrapper = repo_root / "clubb_jax" / "run_jax.py"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_smi = fake_bin / "nvidia-smi"
    fake_smi.write_text(
        "#!/bin/sh\necho 'driver unavailable for test' >&2\nexit 1\n",
        encoding="utf-8",
    )
    fake_smi.chmod(0o755)
    venv = tmp_path / "cuda-venv"
    tools = tmp_path / "tools"
    environment = os.environ | {
        "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        "PYTHON": sys.executable,
        "CLUBB_JAX_VENV": str(venv),
        "CLUBB_JAX_TOOLS_DIR": str(tools),
    }
    environment.pop("CUDA_VISIBLE_DEVICES", None)

    result = subprocess.run(
        [str(wrapper), "--accelerator=cuda13", "--init_env"],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )

    assert result.returncode != 0
    assert "NVIDIA driver is unavailable" in result.stderr
    assert "no CUDA environment was created" in result.stderr
    assert not venv.exists()
    assert not tools.exists()


@pytest.mark.parametrize('device_order', [None, 'FASTEST_FIRST'])
@pytest.mark.parametrize('inherited,flag,expected', [(None, False, 'false'), ('true', False, 'true'), ('false', True, 'true')])
@pytest.mark.parametrize('attached', [False, True])
def test_gpu_wrapper_resolves_physical_index_even_when_pci_order_differs(tmp_path, device_order, inherited, flag, expected, attached):
    repo_root = Path(__file__).resolve().parents[2]
    wrapper = repo_root / "clubb_jax" / "run_jax.py"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_smi = fake_bin / "nvidia-smi"
    fake_smi.write_text(
        "#!/bin/sh\n"
        "echo '0, GPU-test-a, Test GPU A, 610.57.04, 8192, 7.5, 00000000:02:00.0'\n"
        "echo '1, GPU-test-b, Test GPU B, 610.57.04, 16384, 8.6, 00000000:01:00.0'\n",
        encoding="utf-8",
    )
    fake_smi.chmod(0o755)
    environment = os.environ | {
        "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        "PYTHON": sys.executable,
        "CUDA_VISIBLE_DEVICES": "1",
        "CLUBB_JAX_VENV": str(tmp_path / "missing-cuda-venv"),
    }
    environment.pop("CUDA_DEVICE_ORDER", None)
    environment.pop("XLA_PYTHON_CLIENT_PREALLOCATE", None)
    if inherited is not None:
        environment["XLA_PYTHON_CLIENT_PREALLOCATE"] = inherited
    options = ["--accelerator=cuda13"] + (["--xla-prealloc"] if flag else [])
    environment.pop('CLUBB_JAX_REQUESTED_CUDA_VISIBLE_DEVICES', None)
    if device_order:
        environment['CUDA_DEVICE_ORDER'] = device_order

    result = subprocess.run(
        [str(wrapper), *options, "--info=json"],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )

    assert result.returncode == 0
    info = json.loads(result.stdout)
    hardware = info["hardware"]
    assert hardware["cuda_device_order"] == (device_order or "PCI_BUS_ID")
    assert hardware["requested_cuda_visible_devices"] == "1"
    assert hardware["cuda_visible_devices"] == "GPU-test-b"
    assert hardware["selected_gpu"]["index"] == "1"
    assert info["runtime"]["xla_preallocate"] == expected
    assert f"XLA memory preallocation: {expected}" in runtime_info.format_human(info)

    # Observe the environment actually passed into setup, not just the report.
    # Stop at a fake uv so this test never installs packages or runs real GPUs.
    tools_dir = tmp_path / 'tools'
    (tools_dir / 'bin').mkdir(parents=True)
    uv = tools_dir / 'bin' / 'uv'
    uv.write_text('#!/bin/sh\nprintenv CUDA_VISIBLE_DEVICES\nprintenv XLA_PYTHON_CLIENT_PREALLOCATE\nexit 73\n')
    uv.chmod(0o755)
    environment['CLUBB_JAX_TOOLS_DIR'] = str(tools_dir)
    result = subprocess.run(
        [str(wrapper), *options, '--init_env'],
        capture_output=True, text=True, env=environment, timeout=10,
    )
    assert result.returncode == 73
    assert result.stdout.splitlines()[-2:] == ['GPU-test-b', expected]


@pytest.fixture
def gpu_inventory():
    return [
        dict(index='0', uuid='GPU-test-a', name='Test GPU A',
             driver_version='610.10.01', compute_capability='7.0',
             pci_bus_id='0000:02:00.0'),
        dict(index='1', uuid='GPU-test-b', name='Test GPU B',
             driver_version='610.10.01', compute_capability='8.0',
             pci_bus_id='0000:01:00.0'),
    ]


@pytest.mark.parametrize('selector', ['', '-1', '9', 'GPU-missing', 'GPU-test', '0,0', '1,,0'])
def test_invalid_selection_is_rejected(gpu_inventory, selector):
    with pytest.raises(ValueError):
        runtime_info.select_gpus(gpu_inventory, selector)


def test_explicit_selection_preserves_order_and_expands_unique_prefix(gpu_inventory):
    selection = runtime_info.select_gpus(gpu_inventory, 'GPU-test-b,0')
    assert [gpu['index'] for gpu in selection.devices] == ['1', '0']
    assert selection.resolved_devices == 'GPU-test-b,GPU-test-a'


@pytest.mark.parametrize(('selector', 'compatible'), [('0', False), ('1', True), ('1,0', False), ('', False)])
def test_preflight_checks_exposed_gpus_only(monkeypatch, tmp_path, gpu_inventory, selector, compatible):
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', selector)
    monkeypatch.setattr(runtime_info, '_query_nvidia_gpus', lambda: (gpu_inventory, ''))
    requirements = tmp_path / 'requirements.txt'
    requirements.touch()
    info = runtime_info.inspect_runtime('gpu', 'cuda13', requirements, tmp_path / 'venv', '0.11.0')
    assert info['selectable'] is compatible
    assert info['status'] == ('setup_required' if compatible else 'unavailable')


def test_fastest_first_without_selector_does_not_guess_default(monkeypatch, tmp_path, gpu_inventory):
    monkeypatch.delenv('CUDA_VISIBLE_DEVICES', raising=False)
    monkeypatch.setenv('CUDA_DEVICE_ORDER', 'FASTEST_FIRST')
    monkeypatch.setattr(runtime_info, '_query_nvidia_gpus', lambda: (gpu_inventory, ''))
    requirements = tmp_path / 'requirements.txt'
    requirements.touch()
    info = runtime_info.inspect_runtime('gpu', 'cuda13', requirements, tmp_path / 'venv', '0.11.0')
    assert info['hardware']['selected_gpu'] is None
    assert all(gpu['visible'] for gpu in info['hardware']['gpus'])
    assert 'hidden' not in runtime_info.format_human(info)


def test_truncated_nvidia_row_returns_unavailable(monkeypatch):
    monkeypatch.setattr(runtime_info.shutil, 'which', lambda _: '/test/nvidia-smi')
    monkeypatch.setattr(runtime_info.subprocess, 'run', lambda *a, **kw: subprocess.CompletedProcess(
        [], 0, '0,GPU-test,Test GPU,610.10.01,8192,8.0\n', ''
    ))
    gpus, error = runtime_info._query_nvidia_gpus()
    assert not gpus
    assert error


def test_metal_query_accepts_gpu_family_capability(monkeypatch):
    payload = {
        "SPDisplaysDataType": [
            {
                "_name": "Apple M1 Max",
                "sppci_model": "Apple M1 Max",
                "sppci_cores": "32",
                "spdisplays_mtlgpufamilysupport": "spdisplays_metal4",
            }
        ]
    }
    monkeypatch.setattr(runtime_info.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(runtime_info.platform, "machine", lambda: "arm64")
    monkeypatch.setattr(runtime_info.platform, "mac_ver", lambda: ("26.6", ("", "", ""), ""))
    monkeypatch.setattr(
        runtime_info.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess(
            [], 0, json.dumps(payload), ""
        ),
    )

    gpus, error = runtime_info._query_metal_gpus()

    assert error == ""
    assert gpus[0]["name"] == "Apple M1 Max"
    assert gpus[0]["metal_cores"] == 32
    assert gpus[0]["backend"] == "metal"


def test_metal_runtime_uses_plugin_requirements(monkeypatch, tmp_path):
    requirements = tmp_path / "requirements-metal.txt"
    requirements.write_text("jax-metal==0.1.1\n", encoding="utf-8")
    gpu = {
        "index": "0",
        "uuid": "",
        "name": "Apple M1 Max",
        "driver_version": "26.6",
        "memory_mib": None,
        "compute_capability": "",
        "pci_bus_id": "",
        "metal_cores": 32,
        "backend": "metal",
    }
    monkeypatch.setattr(runtime_info, "_query_metal_gpus", lambda: ([gpu], ""))
    monkeypatch.setattr(runtime_info.platform, "mac_ver", lambda: ("26.6", ("", "", ""), ""))

    info = runtime_info.inspect_runtime(
        "gpu", "metal", requirements, tmp_path / "venv", "0.4.34", "3.12"
    )

    assert info["selectable"] is True
    assert info["runtime"]["accelerator"] == "metal"
    assert info["runtime"]["cuda_major"] is None
    assert info["hardware"]["selected_gpu"]["backend"] == "metal"
    assert "Backend: METAL" in runtime_info.format_human(info)
