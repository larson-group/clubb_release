"""Synthetic standalone and input preparation for runner/cohort contracts."""

import json
import sys
from pathlib import Path

from utilities.clubb_case_runner import clubb_case_runner
from input.scripts.create_case_namelist import case_namelist_command


def fake_runner_factory(path, cwd):
    """Return a factory that runs the supplied synthetic standalone in the requested cwd."""
    def factory(args, **kwargs):
        """Build a case runner for one synthetic standalone command."""
        return clubb_case_runner(command=[sys.executable, str(path)], cwd=cwd, **kwargs)
    return factory


def prepare_fake_case(args, output_dir, *, log_file=None, timeout=None):
    """Write JSON containing forwarded CLI arguments and return its synthetic namelist path."""
    namelist = output_dir / f"{args.case_name}.in"
    command = case_namelist_command(args, output_dir)
    namelist.write_text(json.dumps({"argv": [args.case_name, *command[2:-1]]}))
    return namelist


FAKE_STANDALONE = r'''#!/usr/bin/env python3
import json
import os
import sys
import time
from pathlib import Path


def option(name):
    for index, token in enumerate(sys.argv[2:], start=2):
        if token == name:
            return sys.argv[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def run_one(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f" - using executable: {sys.executable}")
    print(" - arguments:", json.dumps(sys.argv[1:]))
    (output_dir / f"{case_name}.in").write_text(
        """&model
    grid_type = 1
    nzmax = 110
    deltaz_nl = 40.0
    zm_init_nl = 0.0
    zm_top_nl = 400.0
    time_initial = 0.0
    time_final = 180.0
    dt_main = 60.0
    /
    """,
        encoding="utf-8",
    )
    (output_dir / f"{case_name}_setup.txt").write_text("setup", encoding="utf-8")
    (output_dir / f"{case_name}_log").write_text("model log", encoding="utf-8")
    started = time.monotonic()
    (output_dir / "invocation.json").write_text(
        json.dumps({"argv": sys.argv[1:], "pid": os.getpid(), "started": started}), encoding="utf-8"
    )

    override = option("-override") or ""
    if "SLEEP=" in override:
        time.sleep(float(override.split("SLEEP=", 1)[1].split(",", 1)[0]))
    else:
        time.sleep(0.15)

    finished = time.monotonic()
    payload = json.loads((output_dir / "invocation.json").read_text(encoding="utf-8"))
    payload["finished"] = finished
    (output_dir / "invocation.json").write_text(json.dumps(payload), encoding="utf-8")
    (Path.cwd() / f"trace_{output_dir.name}.json").write_text(json.dumps(payload), encoding="utf-8")

    if f"FAIL_COLUMNS={columns}" in override:
        print("Synthetic child failure: requested column count is unavailable", file=sys.stderr)
        raise SystemExit(9)

    (output_dir / f"{case_name}.timing").write_text(
        """========================= CLUBB TIMER SUMMARY =========================
    Backend: cpu_time
                                               Timer      Calls   Inclusive(s)   Exclusive(s)     Avg Inc(s)   Excl Pct
    ------------------------------------------------ ---------- -------------- -------------- -------------- ----------
                            advance_clubb_to_end          3       0.300000       0.250000       0.100000      90.00
    ----------------------------------------------------------------------
    Total exclusive timed cpu_time seconds:       0.250000
    =======================================================================
    """,
        encoding="utf-8",
    )

print("CLUBB_WORKER_READY", flush=True)
for command in sys.stdin:
    namelist = Path(command.rstrip("\r\n"))
    sys.argv = [sys.argv[0], *json.loads(namelist.read_text())["argv"]]
    case_name = sys.argv[1]
    columns = int(option("-multicol"))
    run_one(namelist.parent)
    print("CLUBB_WORKER_DONE", flush=True)

'''


def make_fake_standalone(tmp_path: Path) -> Path:
    """Write the synthetic standalone script and return its executable path."""
    path = tmp_path / "fake_standalone.py"
    path.write_text(FAKE_STANDALONE, encoding="utf-8")
    return path
