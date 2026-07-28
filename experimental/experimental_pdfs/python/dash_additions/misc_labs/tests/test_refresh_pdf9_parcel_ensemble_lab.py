"""Command contracts for rebuilding the PDF-9 Notes-lab data."""

from utilities import refresh_pdf9_parcel_ensemble_lab as refresh


def test_refresh_runs_pdf9_cases_then_builds_cache(monkeypatch):
    commands = []
    monkeypatch.setattr(refresh, "_run", lambda command: commands.append(command))

    refresh.refresh()

    assert len(commands) == 3
    for case_name, command in zip(("arm", "bomex"), commands[:2]):
        assert command[2] == "run_scripts/run_scm.py"
        assert command[command.index("-stats") + 1] == str(refresh.STATS_FILE)
        assert command[command.index("-override") + 1] == "iiPDF_type=9"
        assert command[-1] == case_name
    cache_command = commands[2]
    assert cache_command[2] == "utilities/generate_pdf9_parcel_ensemble_diagnostics.py"
    assert cache_command[cache_command.index("--output-dir") + 1] == str(refresh.CACHE_DIR)
    assert f"arm::{refresh.RUNS_DIR / 'arm' / 'arm_stats.nc'}::{refresh.SAM_ROOT / 'arm_3d'}" in cache_command
    assert f"bomex::{refresh.RUNS_DIR / 'bomex' / 'bomex_stats.nc'}::{refresh.SAM_ROOT / 'bomex_3d'}" in cache_command


def test_cache_only_skips_fresh_fortran_runs(monkeypatch):
    commands = []
    monkeypatch.setattr(refresh, "_run", lambda command: commands.append(command))

    refresh.refresh(run_cases=False)

    assert len(commands) == 1
    assert "utilities/generate_pdf9_parcel_ensemble_diagnostics.py" in commands[0]
