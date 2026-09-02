"""Prove the active JAX driver can initialize and step without F2PY."""

from __future__ import annotations

import gc
import importlib.abc
from pathlib import Path
import sys

import jax
import numpy as np
import pytest


class _ClubbPythonBlocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "clubb_python" or fullname.startswith("clubb_python."):
            raise ImportError(f"clubb_python is blocked: {fullname}")
        return None


@pytest.mark.parametrize("case", ["neutral", "arm"])
def test_active_driver_runs_without_clubb_python(case):
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root))
    blocker = _ClubbPythonBlocker()
    sys.meta_path.insert(0, blocker)
    for module_name in tuple(sys.modules):
        if module_name == "clubb_python" or module_name.startswith("clubb_python."):
            del sys.modules[module_name]

    state = None
    try:
        with pytest.raises(ImportError, match="clubb_python is blocked"):
            import clubb_python  # noqa: F401

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
