"""Exercise API-boundary detection and incomplete-input errors on small source trees."""

from pathlib import Path

import pytest

from human_reviewed_tests import check_host_api_imports as checker


@pytest.fixture
def sources(tmp_path):
    """Provide both APIs and private modules whose names differ from their filenames."""
    clubb = tmp_path / "clubb"
    core = clubb / "src/CLUBB_core"
    silhs = clubb / "src/SILHS"
    core.mkdir(parents=True)
    silhs.mkdir()
    (core / "api.F90").write_text("module clubb_api_module\nend module\n")
    (core / "different_filename.F90").write_text(
        "module private_core\nend module\nmodule another_private\nend module\n")
    (silhs / "api.F90").write_text("module silhs_api_module\nend module\n")
    host = tmp_path / "host"
    host.mkdir()
    return clubb, host


@pytest.mark.parametrize("source, line", [
    ("use private_core\nend\n", 1),
    ("program host\n  USE   PRIVATE_CORE, ONLY: x\nend\n", 2),
    ("program host\n  use\tprivate_core\nend\n", 2),
    ("program host\n  use :: private_core\nend\n", 2),
    ("program host\n  use, non_intrinsic :: private_core\nend\n", 2),
    ("program host\n  use &\n! comment\n#ifdef FLAG\n    & private_core\n#endif\nend\n", 2),
    ("program host\n  use private_&\n    &core\nend\n", 2),
    ("program host; use private_core; end program\n", 1),
    ("#if 0\nuse private_core\n#endif\n", 2),
])
def test_host_api_imports_detects_valid_import_forms(sources, source, line):
    """Catch private imports independently of layout, case, continuations and filenames."""
    clubb, host = sources
    path = host / "host.F90"
    path.write_text(source)
    files, modules, violations = checker.check_host_api_imports(clubb, host)
    assert files == 1 and modules == 2
    assert violations == [(path, line, "private_core")]


@pytest.mark.parametrize("source", [
    "use clubb_api_module\nuse silhs_api_module\nuse host_owned_module\nend\n",
    "program host\n! use private_core\nprint *, 'use private_core; ! still text'\nend\n",
    'program host\nprint *, "a doubled ""quote; use private_core"\nend\n',
    'program host\nprint *, "continued &\n&; use private_core"\nend\n',
    "program host\nuse, INTRINSIC :: private_core\nend\n",
])
def test_host_api_imports_ignores_allowed_imports_and_text(sources, source):
    """Allow API/host/intrinsic imports and never treat literal/comment text as code."""
    clubb, host = sources
    (host / "host.f90").write_text(source)
    assert checker.check_host_api_imports(clubb, host)[2] == []


def test_host_api_imports_fixed_form_and_exact_exclusions(sources):
    """Skip vendored source and Git data while checking a similarly named host directory."""
    clubb, host = sources
    vendor = host / "SRC/CLUBB"
    other = host / "other/CLUBB"
    vendor.mkdir(parents=True)
    other.mkdir(parents=True)
    (vendor / "vendor.f90").write_text("use another_private\nend\n")
    path = other / "host.F"
    path.write_text("      program host\nc use another_private\n* use another_private\n"
                    "C use another_private\n      use\n     &private_core\n      end\n")
    (host / ".git").mkdir()
    (host / ".git/hidden.f90").write_text("use another_private\nend\n")
    files, _, violations = checker.check_host_api_imports(clubb, host, ["SRC/CLUBB"])
    assert files == 1
    assert violations == [(path, 5, "private_core")]


@pytest.mark.parametrize("statement", ["      useprivate_core\n",
                                       "      use private_\n     &core\n"])
def test_host_api_imports_fixed_form_blanks(sources, statement):
    """Recognize fixed-form imports with insignificant blanks or split identifiers."""
    clubb, host = sources
    path = host / "host.f"
    path.write_text("      program host\n" + statement + "      end\n")
    assert checker.check_host_api_imports(clubb, host)[2] == [(path, 2, "private_core")]


@pytest.mark.parametrize("problem", ["missing_host", "empty_host", "missing_core", "empty_core",
                                     "missing_api", "bad_exclusion", "exclude_root"])
def test_host_api_imports_incomplete_scan_fails(sources, problem, capsys):
    """Missing/empty source and mistaken exclusions must return input-error status 2."""
    clubb, host = sources
    args = [str(host), "-clubb_source", str(clubb)]
    if problem != "empty_host":
        (host / "host.f90").write_text("use clubb_api_module\nend\n")
    if problem == "missing_host":
        args[0] = str(host / "missing")
    elif problem == "missing_core":
        args[-1] = str(clubb / "missing")
    elif problem == "empty_core":
        for path in (clubb / "src/CLUBB_core").iterdir():
            path.unlink()
    elif problem == "missing_api":
        (clubb / "src/CLUBB_core/api.F90").unlink()
    elif problem == "bad_exclusion":
        args += ["-exclude_dir", "does_not_exist"]
    elif problem == "exclude_root":
        args += ["-exclude_dir", "."]
    assert checker.main(args) == 2
    output = capsys.readouterr()
    assert "Error:" in output.err
    assert "Checked" not in output.out


def test_host_api_imports_legacy_bytes_and_read_failure(sources, monkeypatch, capsys):
    """Find violations after legacy comment bytes and reject an actual read error."""
    clubb, host = sources
    path = host / "host.f90"
    path.write_bytes(b"program host\n! Caf\xe9\nuse private_core\nend\n")
    args = [str(path), "-clubb_source", str(clubb)]
    assert checker.main(args) == 1
    assert f"{path}:3: forbidden import of private_core" in capsys.readouterr().out
    original_open = Path.open

    def unreadable(self, *args, **kwargs):
        """Make only the selected host file unreadable, even under a privileged test user."""
        if self == path:
            raise OSError("simulated unreadable input")
        return original_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "open", unreadable)
    assert checker.main(args) == 2
    assert "simulated unreadable input" in capsys.readouterr().err
