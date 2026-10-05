#!/usr/bin/env python3
"""Initialize and step neutral/ARM with Python API and F2PY imports blocked."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.parse_args()
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()

import gc
import importlib
import importlib.abc
import jax
import numpy as np


class _ClubbPythonBlocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] in {"clubb_python", "clubb_python_api", "clubb_f2py"}:
            raise ImportError(f"clubb_python is blocked: {fullname}")
        return None


def check_active_driver_runs_without_clubb_python(case):
    repo_root = ROOT
    sys.path.insert(0, str(repo_root))
    blocker = _ClubbPythonBlocker()
    sys.meta_path.insert(0, blocker)
    for module_name in tuple(sys.modules):
        if module_name.split(".")[0] in {"clubb_python", "clubb_python_api", "clubb_f2py"}:
            del sys.modules[module_name]

    state = None
    try:
        for name in ("clubb_python", "clubb_python_api", "clubb_f2py"):
            try:
                importlib.import_module(name)
            except ImportError as exc:
                assert "clubb_python is blocked" in str(exc)
            else:
                raise AssertionError(f"Binding import guard failed for {name}")

        from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end
        from clubb_jax.src.clubb_case_initalization import clean_up_clubb, init_clubb_case

        namelist = repo_root / "input" / "case_setups" / f"{case}_model.in"
        state = init_clubb_case(str(namelist))
        advance_clubb_to_end(state, l_stdout=False, max_steps=2)

        thlm = np.asarray(state["thlm"])
        assert thlm.shape == (state["ngrdcol"], state["nzt"])
        assert np.isfinite(thlm).all()
    finally:
        if state is not None:
            clean_up_clubb(state)
        sys.meta_path.remove(blocker)
        if sys.path[0] == str(repo_root):
            sys.path.pop(0)
        jax.clear_caches()
        gc.collect()

if __name__ == "__main__":
    for case in ("neutral", "arm"):
        check_active_driver_runs_without_clubb_python(case)
        print(f"Standalone no-binding validation passed: {case}.")
