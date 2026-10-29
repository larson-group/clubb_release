from types import SimpleNamespace

import pytest

import compile as compiler


@pytest.mark.parametrize("path,expected", [
    ("/usr/bin/gfortran", "gcc"),
    ("/local/bin/x86_64-conda-linux-gnu-gfortran", "gcc"),
    ("/local/bin/aarch64-linux-gnu-gfortran", "gcc"),
    ("/usr/bin/ifx", "intel"),
    ("/usr/bin/nvfortran", "nvhpc"),
])
def test_compiler_executable_maps_to_its_toolchain(path, expected):
    assert compiler.canonical_compiler_from_path(path) == expected


@pytest.mark.parametrize("explicit_fc", [None, "/local/bin/x86_64-conda-linux-gnu-gfortran"])
def test_cmake_respects_explicit_fc_without_overriding_a_toolchain_default(explicit_fc, monkeypatch):
    if explicit_fc:
        monkeypatch.setenv("FC", explicit_fc)
    else:
        monkeypatch.delenv("FC", raising=False)
    monkeypatch.setattr(compiler.shutil, "which", lambda value: explicit_fc if value == explicit_fc else None)
    commands = []
    monkeypatch.setattr(compiler.subprocess, "run", lambda command, **kwargs: commands.append(command))
    args = SimpleNamespace(precision="double", gpu="none", python=False,
                           openmp=False, tuning=False, gptl=False, extra_args=[])
    compiler.configure_cmake(args, "toolchain.cmake", "/install", "Release")
    overrides = [argument for argument in commands[0] if argument.startswith("-DCMAKE_Fortran_COMPILER=")]
    assert overrides == ([f"-DCMAKE_Fortran_COMPILER={explicit_fc}"] if explicit_fc else [])
